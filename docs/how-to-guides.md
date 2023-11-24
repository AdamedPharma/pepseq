This part of the project documentation focuses on a
**problem-oriented** approach. You'll tackle common
tasks that you might have, with the help of the code
provided in this project.

## How To Read Pepseq Format?

You have pepseq format and you need to add them together.
You're in luck! The `pepseq` package can help you
get this done.

Download the code from this GitHub repository and place
the `pepseq/` folder in the same directory as your
Python script:

    your_project/
    │
    ├── pepseq/
    │   ├── __init__.py
    │   └── read.py
    │
    └── your_script.py

Inside of `your_script.py` you can now import the
`read()` function from the `pepseq.read`
module:

    # your_script.py
    from pepseq.read import read

After you've imported the function, you can use it
to add any two numbers that you need to add:

    # your_script.py
    from pepseq.read import read

    print(read('Ac~PSCAG~NH2'))  # OUTPUT: 42.0

You're now able to read pepseq, and you'll
always get a `Peptide` as a result.