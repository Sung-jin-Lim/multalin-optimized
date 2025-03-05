from flask import Flask, request
import subprocess
import os

app = Flask(__name__)

def convert_doc_to_html(filepath):
    """
    Reads an MS Word DOC-format alignment file produced by MultiAlin.
    Assumes that the alignment portion starts after a line beginning with "//".
    Then it converts the special markers:
       - Remove "][" and ")(".
       - Replace "[" with <em class="high"> and "]" with </em>
       - Replace "(" with <em class="low"> and ")" with </em>
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
    
    # Replace markers with HTML tags for coloring.
    text = text.replace("[", '<em class="high">')
    text = text.replace("]", "</em>")
    text = text.replace("(", '<em class="low">')
    text = text.replace(")", "</em>")
    
    # Optionally, wrap the whole thing in a <pre> block for preserving whitespace.
    html = f"<pre>{text}</pre>"
    return html

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
        # Write input to temp_input.fasta
        with open('temp_input.fasta', 'w') as f:
            f.write(raw_sequences)
        
        # Remove previous output DOC files if any, to avoid confusion.
        for file in os.listdir("."):
            if file.startswith("temp_input.doc"):
                os.remove(file)
        
        # Run MultiAlin (ma) with DOC output options.
        # Ensure that your ma binary and .tab files (e.g., blosum62.tab) are accessible.
        # Options: -o:doc forces DOC format; -i:auto autodetects input format; 
        # -k:90.50 sets consensus thresholds; -q quiet; -A aligned.
        subprocess.run([
            './ma',
            '-o:doc',
            'temp_input.fasta'
        ])
        
        # Determine the output file name.
        # Your version of MultiAlin appears to name it "temp_input.doc" (or with numeric suffixes).
        # For simplicity, we try "temp_input.doc" first.
        output_file = "temp_input.doc"
        if not os.path.exists(output_file):
            return """
            <html><body>
            <p>Alignment failed or no output file was created.</p>
            <a href="/">Go Back</a>
            </body></html>
            """
        
        # Convert the DOC file to HTML with color indications.
        alignment_html = convert_doc_to_html(output_file)
        
        return f"""
        <html>
        <head>
          <title>MultiAlin DOC Alignment</title>
          <style>
            pre {{ font-family: monospace; white-space: pre-wrap; }}
            em.high {{ color: red; font-weight: bold; }}
            em.low  {{ color: blue; font-weight: bold; }}
          </style>
        </head>
        <body>
          <h1>Alignment Result (DOC Format with Color Conversion)</h1>
          {alignment_html}
          <br>
          <a href="/">Go Back</a>
        </body>
        </html>
        """
        
    # GET: display input form
    return """
    <html>
    <head><title>Submit Sequences</title></head>
    <body>
      <h1>Submit FASTA Sequences for Alignment</h1>
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
