from flask import Flask, request
import subprocess
import os

app = Flask(__name__)

def convert_doc_to_html(filepath):
    """
    Reads a DOC-format alignment file produced by MultiAlin.
    Assumes the alignment section starts after a line beginning with "//".
    Then it converts special markers:
       - Removes "][" and ")(".
       - Replaces "[" with <em class="high"> and "]" with </em>.
       - Replaces "(" with <em class="low"> and ")" with </em>.
    Returns the resulting HTML string.
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        return f"<p>Error reading DOC file: {e}</p>"
    
    # Find the line that begins with "//"
    start_index = None
    for i, line in enumerate(lines):
        if line.strip().startswith("//"):
            start_index = i + 1
            break
    if start_index is None:
        return "<p>Could not locate alignment section in DOC file.</p>"
    
    # Join the remainder of the file (the alignment block)
    text = "".join(lines[start_index:])
    
    # Remove unwanted adjacent markers
    text = text.replace("][", "")
    text = text.replace(")(", "")
    
    # Replace markers with HTML tags for coloring
    text = text.replace("[", '<em class="high">')
    text = text.replace("]", "</em>")
    text = text.replace("(", '<em class="low">')
    text = text.replace(")", "</em>")
    
    # Wrap the result in a <pre> block
    html = f"<pre>{text}</pre>"
    return html

def parse_msf(filepath):
    """
    A simple parser for MSF formatted files.
    Skips header until a line starting with "//" is found, then processes blocks,
    ignoring lines that start with a digit or "Consensus". Concatenates segments
    for each sequence.
    Returns a list of (name, sequence) tuples.
    """
    sequences = {}
    with open(filepath, 'r') as f:
        lines = f.readlines()
    start_index = 0
    for i, line in enumerate(lines):
        if line.strip().startswith("//"):
            start_index = i + 1
            break
    current_block = []
    for line in lines[start_index:]:
        if not line.strip():
            if current_block:
                for block_line in current_block:
                    parts = block_line.strip().split()
                    if not parts:
                        continue
                    if parts[0].isdigit():
                        continue
                    if parts[0].lower() == "consensus":
                        continue
                    name = parts[0]
                    segment = "".join(parts[1:])
                    sequences.setdefault(name, "")
                    sequences[name] += segment
                current_block = []
        else:
            current_block.append(line)
    if current_block:
        for block_line in current_block:
            parts = block_line.strip().split()
            if not parts:
                continue
            if parts[0].isdigit():
                continue
            if parts[0].lower() == "consensus":
                continue
            name = parts[0]
            segment = "".join(parts[1:])
            sequences.setdefault(name, "")
            sequences[name] += segment
    return [(name, seq) for name, seq in sequences.items()]

def plain_text_alignment(sequences, line_length=60):
    """
    Converts a list of (name, sequence) tuples into plain text.
    Each sequence is printed with its name on the first line, followed by
    the sequence in lines of 'line_length' characters.
    """
    lines = []
    for name, seq in sequences:
        lines.append(name)
        for i in range(0, len(seq), line_length):
            lines.append(seq[i:i+line_length])
        lines.append("")
    return "\n".join(lines)

def parse_editable_alignment(text):
    """
    Parses plain text alignment (snippet style) back into a list of (name, sequence) tuples.
    Expects each sequence block to have the name on the first line followed by sequence lines.
    """
    sequences = []
    blocks = text.strip().split("\n\n")
    for block in blocks:
        lines = block.strip().splitlines()
        if not lines:
            continue
        name = lines[0].strip()
        seq = "".join(line.strip() for line in lines[1:])
        sequences.append((name, seq))
    return sequences

def highlight_columns(seq_chars, high_thresh=0.90, low_thresh=0.50):
    """
    For each column in the alignment (provided as a list of lists of characters),
    applies color coding based on majority rule. Modifies seq_chars in-place.
    """
    import collections
    num_seqs = len(seq_chars)
    if num_seqs == 0:
        return
    align_len = len(seq_chars[0])
    for row in seq_chars:
        if len(row) != align_len:
            return
    for col_idx in range(align_len):
        col = [seq_chars[row_idx][col_idx] for row_idx in range(num_seqs)]
        counter = collections.Counter(col)
        most_common_char, count = counter.most_common(1)[0]
        freq = count / num_seqs
        if most_common_char in ('-', ' '):
            continue
        if freq >= high_thresh:
            for row_idx in range(num_seqs):
                if seq_chars[row_idx][col_idx] == most_common_char:
                    seq_chars[row_idx][col_idx] = f"<em class='high'>{most_common_char}</em>"
        elif freq >= low_thresh:
            for row_idx in range(num_seqs):
                if seq_chars[row_idx][col_idx] == most_common_char:
                    seq_chars[row_idx][col_idx] = f"<em class='low'>{most_common_char}</em>"

def mode3_full_snippet(sequences, chunk_size=60):
    """
    Displays each sequence on its own (name on top, then the sequence in 60-char lines)
    with color-coding applied.
    """
    if not sequences:
        return "<p>No alignment results available.</p>"
    seq_chars = [list(seq) for _, seq in sequences]
    highlight_columns(seq_chars, high_thresh=0.90, low_thresh=0.50)
    align_len = len(seq_chars[0])
    html_out = []
    style_block = """
<style>
pre.seq { font-family: monospace; white-space: pre; }
em.high { color: red; font-weight: bold; }
em.low  { color: blue; font-weight: bold; }
</style>
"""
    html_out.append(style_block)
    html_out.append("<pre class='seq'>")
    for i, (name, _) in enumerate(sequences):
        html_out.append(f"{name}\n")
        modded_seq = seq_chars[i]
        for chunk_start in range(0, align_len, chunk_size):
            chunk_end = min(chunk_start + chunk_size, align_len)
            chunk_slice = "".join(modded_seq[chunk_start:chunk_end])
            html_out.append(chunk_slice + "\n")
        html_out.append("\n")
    html_out.append("</pre>")
    return "".join(html_out)

@app.route('/', methods=['GET', 'POST'])
def index():
    # Main alignment input form
    if request.method == 'POST':
        raw_sequences = request.form.get('sequences', '').strip()
        if not raw_sequences:
            return """
            <html>
            <head>
              <meta name="viewport" content="width=device-width, initial-scale=1">
              <title>No Sequences Provided</title>
              <style>
                body { font-family: Arial, sans-serif; padding: 20px; }
              </style>
            </head>
            <body>
              <p>No sequences provided. Please paste FASTA sequences.</p>
              <a href="/">Go Back</a>
            </body>
            </html>
            """
        with open('temp_input.fasta', 'w') as f:
            f.write(raw_sequences)
        # Remove previous DOC files
        for file in os.listdir("."):
            if file.startswith("temp_input.doc"):
                os.remove(file)
        # Run MultiAlin with DOC output options
        subprocess.run([
            './ma',
            '-o:doc',
            '-r',
            'temp_input.fasta'
        ])
        output_file = "temp_input.doc"
        if not os.path.exists(output_file):
            return """
            <html>
            <head>
              <meta name="viewport" content="width=device-width, initial-scale=1">
              <title>Alignment Failed</title>
              <style>
                body { font-family: Arial, sans-serif; padding: 20px; }
              </style>
            </head>
            <body>
              <p>Alignment failed or no output file was created.</p>
              <a href="/">Go Back</a>
            </body>
            </html>
            """
        alignment_html = convert_doc_to_html(output_file)
        return f"""
        <html>
        <head>
          <meta name="viewport" content="width=device-width, initial-scale=1">
          <title>MultiAlin DOC Alignment</title>
          <style>
            body {{ font-family: Arial, sans-serif; padding: 20px; }}
            pre {{ font-family: monospace; white-space: pre-wrap; }}
            em.high {{ color: red; font-weight: bold; }}
            em.low  {{ color: blue; font-weight: bold; }}
            a {{ display: inline-block; margin-top: 20px; }}
          </style>
        </head>
        <body>
          <h1>Alignment Result (DOC Format)</h1>
          {alignment_html}
          <a href="/">Go Back</a>
          <br>
          <a href="/edit">Edit Alignment</a>
        </body>
        </html>
        """
    return """
    <html>
    <head>
      <meta name="viewport" content="width=device-width, initial-scale=1">
      <title>Submit Sequences</title>
      <style>
        body { font-family: Arial, sans-serif; padding: 20px; }
        textarea { width: 100%; max-width: 600px; height: 200px; font-family: monospace; }
        button { padding: 10px 20px; font-size: 16px; }
      </style>
    </head>
    <body>
      <h1>Submit FASTA Sequences for Alignment</h1>
      <form method="POST">
        <textarea name="sequences" placeholder="Paste FASTA sequences here"></textarea>
        <br><br>
        <button type="submit">Align</button>
      </form>
    </body>
    </html>
    """

@app.route('/edit', methods=['GET', 'POST'])
def edit_alignment():
    if request.method == 'POST':
        edited_text = request.form.get('alignment_text', '').strip()
        if not edited_text:
            return """
            <html>
            <head>
              <meta name="viewport" content="width=device-width, initial-scale=1">
              <title>No Alignment Provided</title>
              <style>
                body { font-family: Arial, sans-serif; padding: 20px; }
              </style>
            </head>
            <body>
              <p>No alignment text provided.</p>
              <a href="/edit">Go Back</a>
            </body>
            </html>
            """
        sequences = parse_editable_alignment(edited_text)
        html_alignment = mode3_full_snippet(sequences, chunk_size=60)
        return f"""
        <html>
        <head>
          <meta name="viewport" content="width=device-width, initial-scale=1">
          <title>Edited Alignment</title>
          <style>
            body {{ font-family: Arial, sans-serif; padding: 20px; }}
            pre {{ font-family: monospace; white-space: pre-wrap; }}
            em.high {{ color: red; font-weight: bold; }}
            em.low  {{ color: blue; font-weight: bold; }}
            a {{ display: inline-block; margin-top: 20px; }}
          </style>
        </head>
        <body>
          <h1>Rechecked Alignment</h1>
          {html_alignment}
          <br>
          <a href="/edit">Edit Again</a>
          <br>
          <a href="/">Go Back to Main</a>
        </body>
        </html>
        """
    else:
        try:
            sequences = parse_msf("temp_input.msf")
        except Exception as e:
            return f"""
            <html>
            <head>
              <meta name="viewport" content="width=device-width, initial-scale=1">
              <title>Error</title>
              <style> body {{ font-family: Arial, sans-serif; padding: 20px; }} </style>
            </head>
            <body>
              <p>Error reading alignment file: {e}</p>
              <a href="/">Go Back</a>
            </body>
            </html>
            """
        plain_text = plain_text_alignment(sequences, line_length=60)
        return f"""
        <html>
        <head>
          <meta name="viewport" content="width=device-width, initial-scale=1">
          <title>Edit Alignment</title>
          <style>
            body {{ font-family: Arial, sans-serif; padding: 20px; }}
            textarea {{ width: 100%; max-width: 600px; height: 300px; font-family: monospace; }}
            button {{ padding: 10px 20px; font-size: 16px; }}
          </style>
        </head>
        <body>
          <h1>Edit Alignment</h1>
          <form method="POST" action="/edit">
            <textarea name="alignment_text">{plain_text}</textarea>
            <br><br>
            <button type="submit">Recheck Alignment</button>
          </form>
          <br>
          <a href="/">Go Back</a>
        </body>
        </html>
        """

if __name__ == '__main__':
    app.run(debug=True)
