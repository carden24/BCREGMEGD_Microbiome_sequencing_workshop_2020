Mohn Lab Best practises for writing bioinformatics scripts
===========================================

### Rationale

Computers are powerful tools for researchers in the life sciences however, scientists are not usually trained in software development. To create reproducible, sound and clear computational research, scientists should follow already established best practises. These practises are oriented toward making research easier to reproduce and scientifically sound.

The following standards are strongly based on the below three articles:

* [Best Practises for Scientific Computing \[DOI: 10.1371/journal.pbio.1001745\]](https://github.com/MohnLab/Lab_Code_Standards_And_Resources/raw/master/Best%20Practices%20for%20Scientific%20Computing.pdf)

* [Ten Simple Rules for Reproducible Computational Research \[DOI: 10.1371/journal.pcbi.1003285\]](https://github.com/MohnLab/Lab_Code_Standards_And_Resources/raw/master/Ten%20Simple%20Rules%20for%20Reproducible%20Computational%20Research.pdf)

* [Ten Simple Rules for the Open Development of Scientific Software \[DOI: 10.1371/journal.pcbi.1002802\]](https://github.com/MohnLab/Lab_Code_Standards_And_Resources/raw/master/Ten%20Simple%20Rules%20for%20the%20Open%20Development%20of%20Scientific%20Software.pdf)

### Guidelines

1. Write programs for people, not for computers. Clarity and readability should be always emphasized
    - Use informative variable names, avoid using names of existing functions (e.g. print) and key words (e.g. import)
    - Variables should be nouns where as functions should be verbs
    - Be consistent with your indentation. Use spaces or tabs but not both
    - Use informative comments that explain what you are doing but do not "state the obvious"
    - Separate names with underscores instead of dots (`variable_1` not `variable.1`). One can also use `CamelCase`.
    - Avoid writing long sections of code, separate the code into multiple blocks (each doing a specific task) 
    - Use meaningful file names and the correct extensions (`.R` not `.r`)
    - Surround binary operators (=, ==, >, <, etc) with whitespaces (`x = y` not `x=y`)
    - Use spaces always after a comma, but never before a comma (`colour = (green, red, blue)` not `colour = (green,red,blue)` or `colour = (green , red , blue)`)
2.	Do not "reinvent the wheel"
    - Check if others have come up with a verified solution to your problem before writing your own scripts
3.	Use standard, well known extensions/packages instead of obscures ones, standard packages are more likely to be tested and maintained
4.	Make incremental changes
    - When writing long code, divide your task in smaller parts that can be performed by custom functions. This is known as [refactoring](http://en.wikipedia.org/wiki/Code_refactoring) your code. Create these functions and test them individually. It will make the whole program/script easier to debug and read
5.	Plan for mistakes
    - Use assertion functions, and internal tests. Do not expect the user to give you the right input
    - Test your scripts with example data, and try to include positive and negative standards. Save the test files
6.	Optimize software only after it works correctly. Work on code readability first, then your algorithm, then optimize.
7.	Document your code's design and purpose, not its mechanics
8.	Backup your data! Store at minimum your raw data and the scripts used to generate the analysis results

### Recommended style guides

For Python, follow [PEP8 standard.](http://legacy.python.org/dev/peps/pep-0008/)
For R use [this standard.](http://stat405.had.co.nz/r-style.html)</p>

### R specific

Save scripts and sessions. If revising old results, open the saved session instead of running the script again, this will save time and avoid changing the results if a random step is involved in your script. 

At the end of your scripts add the following lines to obtain the version of R you are using and the versions of the libraries used:

    R.Version()
    sessionInfo()
    
You can save that information at the end of your final script, or send it directly to a file using:

    writeLines(capture.output(R.Version()), "R_version_info.txt")
    writeLines(capture.output(sessionInfo()), "R_session_info.txt")

### Version Control

Create and use a [Github](https://github.com/) repository to store your working code and keep track its changes. This is helpful for reproducing past results, restoring lost function after drastic code changes, and is essential if you are developing software as a team.

**Note:** Never push your raw data into a repository, only push test data and your scripts!

### Resources

[Software carpentry](http://software-carpentry.org/index.html) has very clear and easy video tutorials on Unix, Python, and version control

### Useful text editors and integrated development environments:

**Notepad++**

This Windows program is very useful to edit codes because it will tell you visually using colors if you have errors in indentation and closing brackets/parenthesis/brace. Download it [here](http://www.notepad-plus-plus.org/)

**TextWrangler**

Rough equivalent of Notepad++ on MAC OSX X. Download it [here](https://itunes.apple.com/ca/app/textwrangler/id404010395?mt=12) 

**Rstudio**

An integrated development environment for developing R scripts and interacting with the R statistical environment. Download it [here](https://www.rstudio.com)
