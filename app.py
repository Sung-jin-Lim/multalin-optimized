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
                # If we already have a sequence buffered, append it
                if name:
                    sequences.append((name, "".join(seq_lines)))
                name = line[1:]  # remove '>'
                seq_lines = []
            else:
                seq_lines.append(line)
        # Append the last sequence if it exists
        if name:
            sequences.append((name, "".join(seq_lines)))

    return sequences

def highlight_alignment(sequences, chunk_size=60, high_thresh=0.90, low_thresh=0.50):
    """
    Takes aligned sequences (list of (name, seq)) and:
      - For each column, finds the majority character frequency
      - If freq >= high_thresh => wrap that char in <em class='high'>
      - Else if freq >= low_thresh => wrap that char in <em class='low'>
      - Else => no <em> (remains black)
    Return an HTML string with line chunking for readability.
    """

    if not sequences:
        return "<p>No alignment results available.</p>"

    num_seqs = len(sequences)
    align_len = len(sequences[0][1])
    # Verify all sequences have the same length
    for name, seq in sequences:
        if len(seq) != align_len:
            return "<p>Error: Inconsistent alignment lengths!</p>"

    # We'll build a 2D array of "HTML-coded characters",
    # so we can fill them column by column.
    # Start by converting each sequence's string into a list of characters:
    seq_chars = [list(seq) for _, seq in sequences]

    # For each column, figure out which character is the majority
    for col_idx in range(align_len):
        # Count frequency of each character in this column
        col = [seq_chars[row_idx][col_idx] for row_idx in range(num_seqs)]
        freq_counter = collections.Counter(col)
        # Most common character and its count
        most_common_char, most_common_count = freq_counter.most_common(1)[0]
        freq = most_common_count / num_seqs  # fraction of sequences that share this char

        # Decide if we do <em class='high'> or <em class='low'>
        if freq >= high_thresh and most_common_char != '-':
            # "High" consensus => red
            # We only wrap the matching characters
            for row_idx in range(num_seqs):
                if seq_chars[row_idx][col_idx] == most_common_char:
                    seq_chars[row_idx][col_idx] = f"<em class='high'>{most_common_char}</em>"
                else:
                    # If a different char or gap, just leave it unwrapped
                    pass

        elif freq >= low_thresh and most_common_char != '-':
            # "Low" consensus => blue
            for row_idx in range(num_seqs):
                if seq_chars[row_idx][col_idx] == most_common_char:
                    seq_chars[row_idx][col_idx] = f"<em class='low'>{most_common_char}</em>"
                # else remain unwrapped

        else:
            # No highlight => remain black text or gap
            pass

    # Build final HTML lines with chunking
    html_out = []
    # Add the style block at top, so we replicate the MultiAlin color style
    style_block = """
<style type='text/css'>
pre.seq  { color: black; background-color: white; }
em.high  { color: red;   background-color: white; font-style: normal; }
em.low   { color: blue;  background-color: white; font-style: normal; }
</style>
"""
    html_out.append(style_block)
    html_out.append("<pre class='seq'>")

    # Now chunk each sequence line by chunk_size
    for i, (name, original_seq) in enumerate(sequences):
        html_out.append(f"{name}\n")  # Print the sequence name on its own line

        # Because we replaced certain chars with <em>..</em>, we need to re-join
        # the modified characters in seq_chars[i].
        modded_seq = seq_chars[i]

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
        # 1) Get sequences from the form
        raw_sequences = request.form.get('sequences', '').strip()
        if not raw_sequences:
            return """
            <html><body>
            <p>No sequences provided. Please paste FASTA sequences.</p>
            <a href="/">Go Back</a>
            </body></html>
            """

        # 2) Write them to a temporary FASTA file
        with open('temp_input.fasta', 'w') as f:
            f.write(raw_sequences)

        # 3) Run Clustal Omega to align (FASTA output)
        subprocess.run([
            'clustalo',
            '-i', 'temp_input.fasta',
            '-o', 'temp_output.aln',
            '--force',
            '--outfmt', 'fa'
        ])

        # 4) Parse the aligned sequences
        aligned_seqs = parse_fasta('temp_output.aln')
        if not aligned_seqs:
            return """
            <html><body>
            <p>Alignment failed or no valid sequences found.</p>
            <a href="/">Go Back</a>
            </body></html>
            """

        # 5) Build color-coded alignment using MultiAlin-like thresholding
        #    90% => red, 50% => blue, else black
        html_alignment = highlight_alignment(aligned_seqs, chunk_size=60,
                                             high_thresh=0.90, low_thresh=0.50)

        # 6) Return the entire HTML page
        return f"""
        <html>
        <head><title>MultiAlin-Style Alignment</title></head>
        <body>
          <h1>MultiAlin-Style Alignment</h1>
          {html_alignment}
          <a href="/">Go Back</a>
        </body>
        </html>
        """

    # If GET, show the input form
    return """
    <html>
    <head><title>Align Sequences</title></head>
    <body>
      <h1>Submit Sequences for MultiAlin-Style Alignment</h1>
      <form method="POST">
        <textarea name="sequences" rows="10" cols="60"
                  placeholder="Paste FASTA sequences here"></textarea>
        <br><br>
        <button type="submit">Align (MultiAlin-Style)</button>
      </form>
    </body>
    </html>
    """

if __name__ == '__main__':
    app.run(debug=True)
