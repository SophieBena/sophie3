
# README for Thermoengine Documentation


## Required installations

- sphinx (1.6.5+, conda-forge)
- numpydoc (conda-forge)
- recommonmark (pip)  
    *Recommonmark allows .md and .rst files in the same project.*

## Structure
The documentation is created using Sphinx and the numpydoc template. It contains mostly rst files but also some legacy md files.

#### Notes
- If you edit a module, you must run `make pyinstall` from the top directory and then `make clean` in the Documentation directory before changes will appear when you run `make html`.

- When you run `make html`, the project throws up a "WARNING: document isn't included in any toctree" for each Markdown file. This is a recommonmark bug/feature.

- Each line in a doc string can be no longer than 80 characters. If you reach 80 characters on a line, simply press Return after the last space, and start on the next line.

- The build configuration file (conf.py) is located at /Documentation/source. Among other things, this file specifies which modules are included in the documentation build.

- Any methods that are inherited from the superclass and that do not have their own documentation will inherit documentation from the superclass. For example, in the Coder module, fdiff methods would inherit the fdiff documentation from the Sympy package.

## Tips

- Make sure not to mix tabs and spaces. If you begin one line with spaces and the others with tabs (for example), you may end up with a blank HTML page or with the inconsistent line omitted.

- If you end up with a blank HTML page, check your indentation. Even one incorrect indent can cause the whole page not to build.

- If you need to update the documentation 


## Equations 
Equations use Mathjax and must be written in LaTeX.

In Markdown files, they must be written like this:

- **Standalone equations** begin with double dollar signs $$ and end with double dollar signs $$

<pre>$$...$$</pre>

- **Inline equations** begin with two backslashes and an open parenthesis `\\(` and end with two backslashes and a close parenthesis `\\)`

``\\(...\\)``

In HTML code, they must be written like this:
- **Standalone equations** are the same as for Markdown.

- **Inline equations** begin with one backslash and an open parenthesis `\(` and end with one backslash and a close parenthesis `\)`

``\(...\)``

In doc strings within modules, they must be written like this:

- **Standalone equations** begin with **:math:\`** and end with **\`**. Make sure that there is not a space after the first tick mark.

- **Inline equations**

For more information about math in modules, see https://documentation.help/Sphinx/math.html




