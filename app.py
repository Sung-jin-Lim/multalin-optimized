from flask import Flask, request
import subprocess
import collections

app = Flask(__name__)

def parse_fasta(filepath):
    """
    Reads a FASTA file and returns a list of (name, sequence) tuples.
    Each sequence is the aligned output from Clustal Omega in FASTA format.
    """
    sequences = []
    name = None
    seq_lines = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name:
                    sequences.append((name, "".join(seq_lines)))
                name = line[1:]  # remove '>'
                seq_lines = []
            else:
                seq_lines.append(line)
        if name:
            sequences.append((name, "".join(seq_lines)))
    return sequences

def highlight_columns(seq_chars, high_thresh=0.90, low_thresh=0.50):
    """
    Given a list of sequences in 'list-of-chars' form (seq_chars),
    apply "majority rule" color-coding:
    - If a char appears in >= high_thresh of sequences => <em class='high'>
    - If >= low_thresh => <em class='low'>
    - Otherwise => black/unwrapped
    - Skip coloring if the majority char is '-'
    Modifies seq_chars in-place.
    """
    import collections

    num_seqs = len(seq_chars)
    if num_seqs == 0:
        return

    align_len = len(seq_chars[0])
    # Verify all sequences have same length
    for row in seq_chars:
        if len(row) != align_len:
            return  # inconsistent length, do nothing

    # For each column, find majority residue and color
    for col_idx in range(align_len):
        col_residues = [seq_chars[row_idx][col_idx] for row_idx in range(num_seqs)]
        counter = collections.Counter(col_residues)
        most_common_char, most_common_count = counter.most_common(1)[0]
        freq = most_common_count / num_seqs

        if most_common_char == '-':
            continue  # skip if the majority is a gap
        if freq >= high_thresh:
            # highlight with <em class='high'>
            for row_idx in range(num_seqs):
                if seq_chars[row_idx][col_idx] == most_common_char:
                    seq_chars[row_idx][col_idx] = f"<em class='high'>{most_common_char}</em>"
        elif freq >= low_thresh:
            # highlight with <em class='low'>
            for row_idx in range(num_seqs):
                if seq_chars[row_idx][col_idx] == most_common_char:
                    seq_chars[row_idx][col_idx] = f"<em class='low'>{most_common_char}</em>"

def mode1_labels_stacked(sequences, chunk_size=60):
    """
    Mode 1: Show "Positions X–Y" before each 60-char chunk,
    and stack sequences below that label.
    """
    if not sequences:
        return "<p>No alignment results available.</p>"

    # Convert each sequence to list-of-chars
    seq_chars = [list(seq) for _, seq in sequences]
    # Color columns in-place
    highlight_columns(seq_chars, high_thresh=0.90, low_thresh=0.50)

    align_len = len(seq_chars[0])
    num_seqs = len(seq_chars)

    html_out = []
    style_block = """
<style>
pre.seq  { color: black; background-color: white; }
em.high  { color: red;   background-color: white; font-style: normal; }
em.low   { color: blue;  background-color: white; font-style: normal; }
</style>
"""
    html_out.append(style_block)
    html_out.append("<pre class='seq'>")

    for chunk_start in range(0, align_len, chunk_size):
        chunk_end = min(chunk_start + chunk_size, align_len)
        # Label
        html_out.append(f"Positions {chunk_start+1}-{chunk_end}\n")
        # Print sequences stacked
        for i, (name, _) in enumerate(sequences):
            # sequence name
            html_out.append(f"{name}\n")
            # chunk slice
            chunk_slice = "".join(seq_chars[i][chunk_start:chunk_end])
            html_out.append(chunk_slice + "\n")
        html_out.append("\n")  # extra spacing

    html_out.append("</pre>")
    return "".join(html_out)

def mode2_no_labels_no_names(sequences, chunk_size=60):
    """
    Mode 2: No "Positions X–Y" labels,
    No sequence names,
    but still chunk at 60 columns and stack sequences together.
    """
    if not sequences:
        return "<p>No alignment results available.</p>"

    seq_chars = [list(seq) for _, seq in sequences]
    highlight_columns(seq_chars, high_thresh=0.90, low_thresh=0.50)

    align_len = len(seq_chars[0])
    num_seqs = len(seq_chars)

    html_out = []
    style_block = """
<style>
pre.seq  { color: black; background-color: white; }
em.high  { color: red;   background-color: white; font-style: normal; }
em.low   { color: blue;  background-color: white; font-style: normal; }
</style>
"""
    html_out.append(style_block)
    html_out.append("<pre class='seq'>")

    for chunk_start in range(0, align_len, chunk_size):
        chunk_end = min(chunk_start + chunk_size, align_len)
        # No label
        for i in range(num_seqs):
            chunk_slice = "".join(seq_chars[i][chunk_start:chunk_end])
            html_out.append(chunk_slice + "\n")
        html_out.append("\n")

    html_out.append("</pre>")
    return "".join(html_out)

def mode3_full_snippet(sequences, chunk_size=60):
    """
    Mode 3: The "snippet style" you provided:
      - One entire sequence at a time
      - Print sequence name
      - Break that sequence at 60 chars
      - Blank line between sequences
      - Column-based color-coding (≥90% => red, ≥50% => blue)
    """

    if not sequences:
        return "<p>No alignment results available.</p>"

    # Convert each sequence to list-of-chars so we can color columns
    seq_chars = [list(seq) for _, seq in sequences]
    highlight_columns(seq_chars, high_thresh=0.90, low_thresh=0.50)

    align_len = len(seq_chars[0])

    html_out = []
    style_block = """
<style type='text/css'>
pre.seq  { color: black; background-color: white; }
em.high  { color: red;   background-color: white; font-style: normal; }
em.low   { color: blue;  background-color: white; font-style: normal; }
</style>
"""
    html_out.append(style_block)
    html_out.append("<pre class='seq'>")

    # Print sequences "one entire sequence at a time" with chunking
    for i, (name, _) in enumerate(sequences):
        html_out.append(f"{name}\n")  # sequence name

        # Re-join the color-coded characters
        modded_seq = seq_chars[i]

        # Break into lines of chunk_size
        for chunk_start in range(0, align_len, chunk_size):
            chunk_end = min(chunk_start + chunk_size, align_len)
            chunk_slice = "".join(modded_seq[chunk_start:chunk_end])
            html_out.append(chunk_slice + "\n")

        html_out.append("\n")  # blank line between sequences

    html_out.append("</pre>")
    return "".join(html_out)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        raw_sequences = request.form.get('sequences', '').strip()
        if not raw_sequences:
            return """
            <html><body>
            <p>No sequences provided. Please paste FASTA sequences.</p>
            <a href="/">Go Back</a>
            </body></html>
            """

        # Run Clustal Omega
        with open('temp_input.fasta', 'w') as f:
            f.write(raw_sequences)
        subprocess.run([
            'clustalo',
            '-i', 'temp_input.fasta',
            '-o', 'temp_output.aln',
            '--force',
            '--outfmt', 'fa'
        ])

        aligned_seqs = parse_fasta('temp_output.aln')
        if not aligned_seqs:
            return """
            <html><body>
            <p>Alignment failed or no valid sequences found.</p>
            <a href="/">Go Back</a>
            </body></html>
            """

        # Build the 3 different HTML displays
        html_mode1 = mode1_labels_stacked(aligned_seqs, chunk_size=60)
        html_mode2 = mode2_no_labels_no_names(aligned_seqs, chunk_size=60)
        html_mode3 = mode3_full_snippet(aligned_seqs, chunk_size=60)

        # Provide a dropdown to switch among them
        final_html = f"""
        <html>
        <head>
          <title>Three Alignment Views</title>
          <script>
          function switchView() {{
            var sel = document.getElementById('viewSelect').value;
            var allViews = ['mode1','mode2','mode3'];
            for(var i=0; i<allViews.length; i++) {{
              document.getElementById(allViews[i]).style.display = 'none';
            }}
            document.getElementById(sel).style.display = '';
          }}
          </script>
        </head>
        <body>
          <h1>Choose Alignment View</h1>
          <p>Select a mode to see how the alignment is displayed:</p>

          <select id="viewSelect" onchange="switchView()">
            <option value="mode1">Mode 1: Labels + Stacked</option>
            <option value="mode2">Mode 2: No Labels, No Names (Stacked)</option>
            <option value="mode3">Mode 3: Snippet-Style (Name + Blank Lines)</option>
          </select>

          <hr>

          <div id="mode1" style="display:;">
            <h2>Mode 1: Labels + Stacked</h2>
            {html_mode1}
          </div>

          <div id="mode2" style="display:none;">
            <h2>Mode 2: No Labels, No Names (Stacked)</h2>
            {html_mode2}
          </div>

          <div id="mode3" style="display:none;">
            <h2>Mode 3: Snippet-Style (One Sequence at a Time)</h2>
            {html_mode3}
          </div>

          <a href="/">Go Back</a>
        </body>
        </html>
        """
        return final_html

    # If GET, show input form
    return """
    <html>
    <head><title>Align Sequences</title></head>
    <body>
      <h1>Submit Sequences for Alignment</h1>
      <form method="POST">
        <textarea name="sequences" rows="10" cols="60"
                  placeholder="Paste FASTA sequences here"></textarea>
        <br><br>
        <button type="submit">Align</button>
      </form>
    </body>
    </html>
    """

if __name__ == '__main__':
    app.run(debug=True)
