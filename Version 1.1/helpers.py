import os
import requests
import urllib.parse

from flask import redirect, render_template, request, session
from functools import wraps


def transform_view(file_name):
    file = request.files['data_file']
    if not file:
        return "No file"

    file_contents = list(file.stream.read().decode("utf-8").splitlines())
    reads_freq = collections.Counter(file_contents)
    seqs = list(reads_freq.keys())
    counts = list(reads_freq.values())

    round = open(file_name).readlines()

    quality = round[3::4]
    seqs = round[1::4]

    return """
            <html>
                <body>
                   Sequences: {seqs}
                   Counts: {counts}
                </body>
            </html>
        """.format(counts = counts, seqs = seqs)

def error_page(message, code=777):
    """Render message as an apology to user."""
    def escape(s):
        """
        Escape special characters.

        https://github.com/jacebrowning/memegen#special-characters
        """
        for old, new in [("-", "--"), (" ", " "), ("_", "__"), ("?", "~q"),
                         ("%", "~p"), ("#", "~h"), ("/", "~s"), ("\"", "''")]:
            s = s.replace(old, new)
        return s
    return render_template("error_page.html", err=escape(message)), code
