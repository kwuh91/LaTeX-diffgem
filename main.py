import sympy   as sp
import numpy   as np
import extra

from LaTeXpy import LaTeX


# ⁱᵢʲⱼ

if __name__ == '__main__':
    # packages
    packages = { # from : what
        'inputenc'  : ['utf8'],  
        'babel'     : ['russian'],
        'amsmath'   : [],
        'esdiff'    : ['thinc'],
        'mathtools' : [],
        'fontenc'   : ['T1'],
    }

    # create file for output and LaTeX instance
    with open("latex.tex", "w", encoding="utf-8") as f, LaTeX(f, packages=packages) as doc:
        
        m   = LaTeX.m   # math mode
        bf  = LaTeX.bf  # bold font   
        rf  = LaTeX.rf  # roman font  
        sup = LaTeX.sup # superscript 
        sub = LaTeX.sub # subscript   
        eq  = LaTeX.eq  # equation
        lb  = LaTeX.lb  # line break 

        bs  = '\\' # backslash

        X, Y, Z = sp.symbols(f'X¹ X² X³')

        # связь с цилиндрическими координатами
        cyl_coord = {
            'x¹' : X * sp.cos(Y),
            'x²' : X * sp.sin(Y),
            'x³' : Z,
        }

        # символ кронекера
        δ = [[1, 0, 0],
             [0, 1, 0],
             [0, 0, 1]]

        Tⁱʲ = [[0,   X-Y, 0],
               [0,   0,   0],
               [2*Z, 0,   0]]

        doc.section('№3')
        doc.subsection('Условие:')
        doc.write(f'Задано тензорное поле {m(rf("T"))}({m("X^i")}), где {m("X^i")} - цилиндрические координаты. Найти:')
        doc.write(f'{m(rf("1)"))} ковариантные, контравариантные компоненты этого поля в базисах {m(rf("r_i"))}, {m(rf("r^i"))},', end='')
        doc.write(f'где {m(rf("r_i"))} — ортогональный локальный базис цилиндрической системы координат;')
        doc.write(f'{m(rf("2)"))} ковариантные производные компонент тензорного поля в базисах {m(rf("r_i"))}, {m(rf("r^i"))}')

        doc.subsection('Исходные данные:')
        doc.write(f'Поле задано в виде: {m(rf("T"))} = {m("T^ij")}{m(rf("e_i"))}$\\otimes${m(rf("e_j"))};')

        doc.write(f'\n{eq(f"T^{{ij}}={LaTeX.plainMatrix3D(Tⁱʲ)}")}', end=f'\n')

        doc.subsection('Решение:')
        doc.write(f'{m(rf("r_i"))} - ортогональный локальный базис цилиндрической системы координат:')

        # cyl coordinates
        first  = f'x^1 = X^1 \\cdot \\cos(X^2)'
        second = f'x^2 = X^1 \\cdot \\sin(X^2)'
        third  = f'x^3 = X^3'

        doc.write(f'{eq(LaTeX.system(first, second, third))}', end=f'{lb}')

        doc.write(f'{m(rf("1)"))} Для того, чтобы найти компоненты тензорного поля в базисах {m(rf("r_i"))}, {m(rf("r^i"))},', end='')
        doc.write(f'найдем сначала метрическую матрицу для цилиндрических координат {m("X^i")} и локальные векторы базиса.') 

        doc.write(f'{m(rf("1.1)"))} Найдем якобиеву матрицу для криволинейных координат {m("X^i")}:')

        doc.write(eq(f'Q^i_k = {bs}frac{{{bs}partial x^i}}{{{bs}partial X^k}}'), end='')

        # якобиева матрица цилиндрической системы координат
        QⁱₖSolutions = [[0 for _ in range(3)] for _ in range(3)]
        Qⁱₖ          = [[0 for _ in range(3)] for _ in range(3)]

        # Вычислим её компоненты
        for i in range(1, 3 + 1):
            for k in range(1, 3 + 1):

                res = sp.diff(cyl_coord[f'x{extra.get_super(i)}'], f'X{extra.get_super(k)}')
                resStr = LaTeX.cleanUp(res)
                
                QⁱₖSolutions[i - 1][k - 1] = f'Q^{i}_{k} = {bs}frac{{{bs}partial x^{i}}}{{{bs}partial X^{k}}} = {resStr}'
                Qⁱₖ[i - 1][k - 1]          = res

        doc.write(f'{LaTeX.alignat(QⁱₖSolutions, quad = True)}')

        doc.write(f'Тогда якобиева матрица цилиндрической системы координат имеет вид:')
        doc.write(f'\n{eq(f"Q^i_k={LaTeX.plainMatrix3D(Qⁱₖ)}")}', end=f'{lb}')

        doc.write(f'{m(rf("1.2)"))} Найдем локальные векторы базиса для цилиндрических координат {m("X^i")}:')
        doc.write(eq(f'{rf("r_k")} = {bs}frac{{{bs}partial x^i}}{{{bs}partial X^k}}\\overline{{{rf("e_k")}}} = Q^i_k\\overline{{{rf("e_k")}}}'))
        doc.write(f'Матрица {m("Q^i_k")} была найдена на предыдущем шаге.')

        doc.write(f'{m(rf("1.3)"))} Найдем метрическую матрицу для цилиндрических координат {m("X^i")}:')
        doc.write(eq(f'{rf("g_{{ij}}")} = {rf("r_i")}\\cdot{rf("r_j")} = Q^s_iQ^p_j\\delta_{{sp}}'))

        doc.write(f'Найдем компоненты метрической матрицы для цилиндрической системы координат:')

        # метрическая матрица для цилиндрической системы координат
        gSolutions = [0 for _ in range(3*3)]
        g          = [[0 for _ in range(3)] for _ in range(3)]

        count = 0
        
        # Вычислим её компоненты
        for i in range(1, 3 + 1):
            for j in range(1, 3 + 1):
                
                res = 0
                resStr = rf(f"g_{{{i}{j}}}")

                resStr += f' = Q^s_{i} \\cdot Q^p_{j} * \\delta_{{sp}} = ' 
                
                for s in range(1, 3 + 1):
                    for p in range(1, 3 + 1):

                        mid_res = Qⁱₖ[s - 1][i - 1] * Qⁱₖ[p - 1][j - 1] * δ[s - 1][p - 1]
                        res += mid_res
                
                res     = sp.simplify(res)
                resStr += LaTeX.cleanUp(res)

                g[i - 1][j - 1]   = res
                gSolutions[count] = resStr
                count += 1                

        doc.write(f'{LaTeX.alignat(gSolutions, quad = True)}')

        doc.write(f'Запишем полученную метрическую матрицу для цилиндрической системы координат:')
        doc.write(f'{eq(f"g_{{ij}}={LaTeX.plainMatrix3D(g)}")}', end=f'{lb}')

        doc.write(f'Найдем обратную метрическую матрицу:')
        gInverse = extra.getMatrixInverse(g)
        doc.write(f'{eq(f"g^{{ij}}={LaTeX.plainMatrix3D(gInverse)}")}', end=f'{lb}')

        doc.write(f'{m(rf("1.4)"))} Найдем векторы взаимного локального базиса для цилиндрических координат {m("X^i")}:')

        doc.write(eq(f'{rf("r^{i}")} = {rf("g^{{ij}}")}{rf("r_j")} = {rf("g^{{ij}}")}Q^m_j\\overline{{{rf("e_m")}}} = Q^{{im}}\\overline{{{rf("e_m")}}}:'))

        QⁱᵐSolutions = [0 for _ in range(3*3)]
        Qⁱᵐ          = [[0 for _ in range(3)] for _ in range(3)]

        count = 0
        
        # Вычислим её компоненты
        for i in range(1, 3 + 1):
            for m in range(1, 3 + 1):
                
                res = 0
                resStr = ''
                midResArr = []

                resStr += f'Q^{{{i}{m}}} = '
                resStr += rf(f'g^{{{i}j}}')
                resStr += f'Q^{m}_j = '

                for j in range(1, 3 + 1):

                    midRes = gInverse[i - 1][j - 1] * Qⁱₖ[m - 1][j - 1]
                    midResArr.append(midRes)

                    res += midRes
                    resStr += f'g^{{{i}{j}}} \\cdot Q^{m}_{j}'
                    resStr += f'+ ' if j < 3 - 1 else f'= '
                
                for it in range(len(midResArr)):
                    resStr += LaTeX.cleanUp(midResArr[it])
                    resStr += f'+ ' if it < len(midResArr) - 1 else f'= '

                res     = sp.simplify(res)
                resStr += LaTeX.cleanUp(res)

                Qⁱᵐ[i - 1][m - 1]   = res
                QⁱᵐSolutions[count] = resStr
                count += 1    
        
        doc.write(f'{LaTeX.alignat(QⁱᵐSolutions, quad = True)}')

        doc.write(f'Запишем полученную матрицу:')
        doc.write(f'{eq(f"Q^{{im}}={LaTeX.plainMatrix3D(Qⁱᵐ)}")}', end=f'{lb}')

        m = LaTeX.m 
        doc.write(f'{m(rf("1.5)"))}Найдем теперь компоненты тензорного поля.')
