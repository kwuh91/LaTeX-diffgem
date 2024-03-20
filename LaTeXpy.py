from typing import List
import re

class LaTeX:
    # Latex constructor
    def __init__(self, file, packages    = None, 
                             title       = None, 
                             author      = None, 
                             date        = None):
        # assign file fields
        self.file = file

        # initialize doc
        self.file.write("\\documentclass{article}\n")

        # import packages
        if packages is not None:

            for package_name, value in packages.items():
                size = len(value) # amount of explicitly chosen names from a package

                if size > 0:
                    self.file.write(f"\\usepackage[")

                    for i in range(len(value)):
                        self.file.write(f"{value[i]}")
                        self.file.write(", ") if i < len(value) - 1 else self.file.write(f"]{{{package_name}}}\n")
                else:
                    self.file.write(f"\\usepackage{{{package_name}}}\n")
        
        # line skip
        self.blank()

        # write title
        if title is not None:
            self.file.write(f"\\title{{{title}}}\n")

        # write author
        if author is not None:
            self.file.write(f"\\author{{{author}}}\n")

        # write date
        if title        is not None and \
           author       is not None and \
           date         is not None:
           self.file.write(f"\\date{{{date}}}\n")

        # add line skip if needed
        if title  is not None and \
           author is not None:
            self.blank()

        self.file.write("\\begin{document}\n\n")

        self.file.write("\\vspace*{-0cm} % Adjust the value as needed to remove the white space\n\n")

        # create title if needed
        if title  is not None and \
           author is not None:
            self.file.write("\\maketitle\n")

    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_value, tb):
        self.file.write("\\vspace*{-0cm} % Adjust the value as needed to remove the white space\n\n")
        self.file.write("\\end{document}\n")
        self.file.close()

    m    = staticmethod(lambda a : f'${a}$')            # math mode
    bf   = staticmethod(lambda a : f'\mathbf{{{a}}}')   # bold font
    rf   = staticmethod(lambda a : f'\mathrm{{{a}}}')   # roman font
    sup  = staticmethod(lambda a : f'^{a}')             # superscript
    sub  = staticmethod(lambda a : f'_{a}')             # subscript
    eq   = staticmethod(lambda a : f'\\[\n{a}\n\\]')    # equation
    lb   = '\\\\'                                       # line break 

    indent  = '  '                                      # indent with the size of 2 spaces
    cleanUp = staticmethod(lambda a : LaTeX.replace_double_parentheses(
                                      LaTeX.replace_substrings(str(a).replace('X¹','(X^1)')  \
                                                                     .replace('X²','(X^2)')  \
                                                                     .replace('X³','(X^3)')  \
                                                                     .replace('cos','\\cos') \
                                                                     .replace('sin','\\sin') \
                                                                     .replace('**', '^')     \
                                                                     .replace('*', '\\cdot '))))

    @staticmethod
    def plainMatrix3D(matrix: List[List[any]]):        # matrix 
        res = f'\\begin{{pmatrix}}\n'                  # Parentheses;
        for i in range(len(matrix)):                   # round brackets 
            res += f'\t'
            for j in range(len(matrix[i])):

                res += LaTeX.cleanUp(matrix[i][j])
                
                if j < len(matrix[i]) - 1:
                    res += f' & '
                elif i < len(matrix[i]) - 1:
                    res += f'{LaTeX.lb}\n'
                    
        res += f'\n\end{{pmatrix}}'
        return res

    @staticmethod                                       # cases
    def system(*args: List[any]):
        res = f'\\begin{{cases}}\n'
        for i in range(len(args)):
            res += f'{LaTeX.indent}{args[i]}'
            res += f'\\\\\n' if i < len(args) - 1 else f'\n'
        res += f'\\end{{cases}}'
        return res

    @staticmethod # [[any for _ in range (len(arr[_]))] for _ in range (any)]
    def alignat(args: List[List[any]], quad = False): # alignment
        if len(args) > 0:
            res = ''
            rowSize = 0
            twoDim = False

            if isinstance(args, list)      and \
               all((isinstance(row, list)) for row in args):
                rowSize = len(args[0])
                twoDim = True

            elif isinstance(args, list)    and \
                 all((isinstance(row, str)) for row in args):
                rowSize = 1

            res = f'\\begin{{alignat*}}{{{rowSize}}}\n'
            quad = '\\quad' if quad else ''

            if twoDim: # 2D args
                for i in range(len(args)):
                    res += LaTeX.indent + '& '
                    for j in range(rowSize):
                        res += f'{args[i][j]} '

                        if j < rowSize - 1:
                            res += f'{quad} &&'
                        elif i < len(args) - 1:
                            res += LaTeX.lb + '\n'
                        else:
                            res += f'\n'
            else: # 1D args
                for i in range(len(args)):
                    res += LaTeX.indent + f'& {args[i]} '

                    if i < len(args) - 1:
                        res += LaTeX.lb + '\n'
                    else:
                        res += f'\n'

            res += f'\\end{{alignat*}}'
        else:
            raise Exception('Invalid parameter\\s')
    
        return res

    @staticmethod # turn "(...)^(-...)" into "//frac{1}{(...)^(...)}"
    def replace_substrings(input_string):
        # Define the regular expression pattern to match "(...)^(-...)"
        pattern = r"\((.*?)\)\^\((-.*?)\)"

        def repl(match):
            base = match.group(1)
            exponent = match.group(2)
            # Check if exponent is negative
            if exponent.startswith("-"):
                return f"\\frac{{1}}{{({base})^{exponent[1:]}}}"
            else:
                return f"\\frac{{1}}{{({base})^{exponent}}}"

        # Use re.sub to replace substrings matching the pattern
        result = re.sub(pattern, repl, input_string)

        return result

    @staticmethod # This will replace all occurrences of '((...))' with '(...)' in the given text.
    def replace_double_parentheses(text):
        # return re.sub(r'\(\((.*?)\)\)', r'(\1)', text)
        return text

    # skip line
    def blank(self):
        self.file.write('\n')

    # create a section
    def section(self, text, numbering=False):
        if numbering:
            self.file.write(f"\\section{{{text}}}\n")
        else:
            self.file.write(f"\\section*{{{text}}}\n")

    # create a section
    def subsection(self, text, numbering=False):
        if numbering:
            self.file.write(f"\\subsection{{{text}}}\n")
        else:
            self.file.write(f"\\subsection*{{{text}}}\n")

    # append text to Latex document
    def write(self, text, end=None):
        end = '\\\\\n' if end is None else end + '\n'
        self.file.write(text + end)

    # Latex destructor
    def __del__(self):
        pass
