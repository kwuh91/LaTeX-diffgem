import sympy   as sp
import numpy   as np
import extra

from typing import List, Any
from LaTeXpy import LaTeX

# ⁱᵢʲⱼ

# Расчет ковариантных производных тензорного поля
def covariant(T: List[List[Any]], Gamma: List[List[List[Any]]], var: int):
        
        m   = LaTeX.m   # math mode 
        rf  = LaTeX.rf  # roman font  
        eq  = LaTeX.eq  # equation
        lb  = LaTeX.lb  # line break 

        bs  = '\\' # backslash

        match var:
            case 1: # nabla_k T^{ij}

                doc.write(m(rf(f'2.{var}) ')) + f'Вычислим ковариантную производную контравариантных компонент поля по формуле:')
                doc.write(eq(f'{bs}nabla_kT^{{ij}} = {bs}frac{{{bs}partial T^{{ij}}}}{{{bs}partial X^k}} + T^{{mj}}{bs}Gamma^i_{{mk}} + T^{{im}}{bs}Gamma^j_{{mk}};'))

            case 2: # nabla_k T_{ij}

                doc.write(m(rf(f'2.{var}) ')) + f'Вычислим ковариантную производную ковариантных компонент поля по формуле:')
                doc.write(eq(f'{bs}nabla_kT_{{ij}} = {bs}frac{{{bs}partial T_{{ij}}}}{{{bs}partial X^k}} - T_{{mj}}{bs}Gamma^m_{{ik}} - T_{{im}}{bs}Gamma^m_{{jk}};'))

            case 3: # nabla_k T^i_j

                doc.write(m(rf(f'2.{var}) ')) + f'Вычислим ковариантную производную смешанных компонент поля по формуле:')
                doc.write(eq(f'{bs}nabla_kT^i_j = {bs}frac{{{bs}partial T^i_j}}{{{bs}partial X^k}} + T^m_j{bs}Gamma^i_{{mk}} - T^i_m{bs}Gamma^m_{{jk}};'))

            case 4: # nabla_k T^j_i

                doc.write(m(rf(f'2.{var}) ')) + f'Вычислим ковариантную производную смешанных компонент поля по формуле:')
                doc.write(eq(f'{bs}nabla_kT^j_i = {bs}frac{{{bs}partial T^j_i}}{{{bs}partial X^k}} - T^j_m{bs}Gamma^m_{{ik}} + T^m_i{bs}Gamma^j_{{mk}};'))

        # контейнер для ковариантной производной контравариантных компонент
        nablaSolutions = [[0 for _ in range(3 * 3)] for _ in range(3)] 
        nabla          = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)]

        for k in range(1, 3 + 1):
            count = 0
            for i in range(1, 3 + 1):
                for j in range(1, 3 + 1):
                    
                    res    = 0
                    resStr = ''

                    match var:
                        case 1: # nabla_k T^{ij}

                            resStr += f'\\nabla_{k}T^{{{i}{j}}} = '
                            resStr += f'{bs}frac{{{bs}partial T^{{{i}{j}}}}}{{{bs}partial X^{k}}} + '
                            resStr += f'T^{{m{j}}}{bs}Gamma^{i}_{{m{k}}} + '
                            resStr += f'T^{{{i}m}}{bs}Gamma^{j}_{{m{k}}} = '

                            res += sp.diff(T[i - 1][j - 1], f'X{extra.get_super(k)}')
                            for m in range(1, 3 + 1):
                                res = res + \
                                      T[m - 1][j - 1] * Gamma[i - 1][m - 1][k - 1] + \
                                      T[i - 1][m - 1] * Gamma[j - 1][m - 1][k - 1]

                        case 2: # nabla_k T_{ij}
                            
                            resStr += f'\\nabla_{k}T_{{{i}{j}}} = '
                            resStr += f'{bs}frac{{{bs}partial T_{{{i}{j}}}}}{{{bs}partial X^{k}}} - '
                            resStr += f'T_{{m{j}}}{bs}Gamma^m_{{{i}{k}}} - '
                            resStr += f'T_{{{i}m}}{bs}Gamma^m_{{{j}{k}}} = '

                            res += sp.diff(T[i - 1][j - 1], f'X{extra.get_super(k)}')
                            for m in range(1, 3 + 1):
                                res = res - \
                                      T[m - 1][j - 1] * Gamma[m - 1][i - 1][k - 1] - \
                                      T[i - 1][m - 1] * Gamma[m - 1][j - 1][k - 1]

                        case 3: # nabla_k T^i_j
                            
                            resStr += f'\\nabla_{k}T^{i}_{j} = '
                            resStr += f'{bs}frac{{{bs}partial T^{i}_{j}}}{{{bs}partial X^{k}}} + '
                            resStr += f'T^m_{j}{bs}Gamma^{i}_{{m{k}}} - '
                            resStr += f'T^{i}_m{bs}Gamma^m_{{{j}{k}}} = '

                            res += sp.diff(T[i - 1][j - 1], f'X{extra.get_super(k)}')
                            for m in range(1, 3 + 1):
                                res = res + \
                                      T[m - 1][j - 1] * Gamma[i - 1][m - 1][k - 1] - \
                                      T[i - 1][m - 1] * Gamma[m - 1][j - 1][k - 1]

                        case 4: # nabla_k T^j_i
                            
                            resStr += f'\\nabla_{k}T^{j}_{i} = '
                            resStr += f'{bs}frac{{{bs}partial T^{j}_{i}}}{{{bs}partial X^{k}}} - '
                            resStr += f'T^{j}_m{bs}Gamma^m_{{{i}{k}}} + '
                            resStr += f'T^m_{i}{bs}Gamma^{j}_{{m{k}}} = '

                            res += sp.diff(T[i - 1][j - 1], f'X{extra.get_super(k)}')
                            for m in range(1, 3 + 1):
                                res = res - \
                                      T[m - 1][j - 1] * Gamma[m - 1][i - 1][k - 1] + \
                                      T[i - 1][m - 1] * Gamma[j - 1][m - 1][k - 1]

                    res     = sp.simplify(res)
                    resStr += LaTeX.cleanUp(res)
                    
                    nabla[k - 1][i - 1][j - 1]   = res
                    nablaSolutions[k - 1][count] = resStr
                    count += 1

        doc.write(f'При k = 1:')
        doc.write(f'{LaTeX.alignat(nablaSolutions[1 - 1], quad = True)}')

        doc.write(f'При k = 2:')
        doc.write(f'{LaTeX.alignat(nablaSolutions[2 - 1], quad = True)}')

        doc.write(f'При k = 3:')
        doc.write(f'{LaTeX.alignat(nablaSolutions[3 - 1], quad = True)}')

        doc.write(f'Запишем результат:')

        match var:
            case 1: # nabla_k T^{ij}

                doc.write(f'{eq(f"{bs}nabla_1T^{{ij}} = {LaTeX.plainMatrix3D(nabla[1 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla_2T^{{ij}} = {LaTeX.plainMatrix3D(nabla[2 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla_3T^{{ij}} = {LaTeX.plainMatrix3D(nabla[3 - 1])}")}', end=f'{lb}')

            case 2: # nabla_k T_{ij}

                doc.write(f'{eq(f"{bs}nabla_1T_{{ij}} = {LaTeX.plainMatrix3D(nabla[1 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla_2T_{{ij}} = {LaTeX.plainMatrix3D(nabla[2 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla_3T_{{ij}} = {LaTeX.plainMatrix3D(nabla[3 - 1])}")}', end=f'{lb}')

            case 3: # nabla_k T^i_j

                doc.write(f'{eq(f"{bs}nabla_1T^i_j = {LaTeX.plainMatrix3D(nabla[1 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla_2T^i_j = {LaTeX.plainMatrix3D(nabla[2 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla_3T^i_j = {LaTeX.plainMatrix3D(nabla[3 - 1])}")}', end=f'{lb}')

            case 4: # nabla_k T^j_i
                doc.write(f'{eq(f"{bs}nabla_1T^j_i = {LaTeX.plainMatrix3D(nabla[1 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla_2T^j_i = {LaTeX.plainMatrix3D(nabla[2 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla_3T^j_i = {LaTeX.plainMatrix3D(nabla[3 - 1])}")}', end=f'{lb}')

        return [nablaSolutions, nabla]


# Расчет контравариантных производных тензорного поля
def contravariant(gI: List[List[Any]], nabla: List[List[List[Any]]], var: int):
        
        m   = LaTeX.m   # math mode 
        rf  = LaTeX.rf  # roman font  
        eq  = LaTeX.eq  # equation
        lb  = LaTeX.lb  # line break 

        bs  = '\\' # backslash

        match var:
            case 1: # nabla^m T^{ij}

                doc.write(m(rf(f'2.{var + 4}) ')) + f'Вычислим контравариантную производную от контравариантных компонент по формуле:')
                doc.write(eq(f'{bs}nabla^mT^{{ij}} = g^{{mk}}{bs}nabla_kT^{{ij}};'))

            case 2: # nabla^m T_{ij}

                doc.write(m(rf(f'2.{var + 4}) ')) + f'Вычислим контравариантную производную от ковариантных компонент по формуле:')
                doc.write(eq(f'{bs}nabla^mT_{{ij}} = g^{{mk}}{bs}nabla_kT_{{ij}};'))

            case 3: # nabla^m T^i_j

                doc.write(m(rf(f'2.{var + 4}) ')) + f'Аналогично вычислим контравариантную производную от смешанных компонент тензора по формуле:')
                doc.write(eq(f'{bs}nabla^mT^i_j = g^{{mk}}{bs}nabla_kT^i_j;'))

            case 4: # nabla_k T^j_i

                doc.write(m(rf(f'2.{var + 4}) ')) + f'Аналогично вычислим контравариантную производную от смешанных компонент тензора по формуле:')
                doc.write(eq(f'{bs}nabla^mT^j_i = g^{{mk}}{bs}nabla_kT^j_i;'))

        # контейнер для ковариантной производной контравариантных компонент
        nablaSolutions2 = [[0 for _ in range(3 * 3)] for _ in range(3)] 
        nabla2          = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)]

        for m in range(1, 3 + 1):
            count = 0
            for i in range(1, 3 + 1):
                for j in range(1, 3 + 1):
                    
                    res    = 0
                    resStr = ''

                    match var:
                        case 1: # nabla^m T^{ij}

                            resStr += f'\\nabla^{m}T^{{{i}{j}}} = '
                            resStr += f'g^{{{m}k}}{bs}nabla_kT^{{{i}{j}}} = '

                        case 2: # nabla^m T_{ij}
                            
                            resStr += f'\\nabla^{m}T_{{{i}{j}}} = '
                            resStr += f'g^{{{m}k}}{bs}nabla_kT_{{{i}{j}}} = '

                        case 3: # nabla^m T^i_j
                            
                            resStr += f'\\nabla^{m}T^{i}_{j} = '
                            resStr += f'g^{{{m}k}}{bs}nabla_kT^{i}_{j} = '

                        case 4: # nabla^m T^j_i
                            
                            resStr += f'\\nabla^{m}T^{j}_{i} = '
                            resStr += f'g^{{{m}k}}{bs}nabla_kT^{j}_{i} = '

                    for k in range(1, 3 + 1):
                        res = res + \
                                gI[m - 1][k - 1] * nabla[k - 1][i - 1][j - 1]


                    res     = sp.simplify(res)
                    resStr += LaTeX.cleanUp(res)
                    
                    nabla2[m - 1][i - 1][j - 1]   = res
                    nablaSolutions2[m - 1][count] = resStr
                    count += 1

        doc.write(f'При k = 1:')
        doc.write(f'{LaTeX.alignat(nablaSolutions2[1 - 1], quad = True)}')

        doc.write(f'При k = 2:')
        doc.write(f'{LaTeX.alignat(nablaSolutions2[2 - 1], quad = True)}')

        doc.write(f'При k = 3:')
        doc.write(f'{LaTeX.alignat(nablaSolutions2[3 - 1], quad = True)}')

        doc.write(f'Запишем результат:')

        match var:
            case 1: # nabla^m T^{ij}

                doc.write(f'{eq(f"{bs}nabla^1T^{{ij}} = {LaTeX.plainMatrix3D(nabla2[1 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla^2T^{{ij}} = {LaTeX.plainMatrix3D(nabla2[2 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla^3T^{{ij}} = {LaTeX.plainMatrix3D(nabla2[3 - 1])}")}', end=f'{lb}')

            case 2: # nabla^m T_{ij}

                doc.write(f'{eq(f"{bs}nabla^1T_{{ij}} = {LaTeX.plainMatrix3D(nabla2[1 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla^2T_{{ij}} = {LaTeX.plainMatrix3D(nabla2[2 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla^3T_{{ij}} = {LaTeX.plainMatrix3D(nabla2[3 - 1])}")}', end=f'{lb}')

            case 3: # nabla^m T^i_j

                doc.write(f'{eq(f"{bs}nabla^1T^i_j = {LaTeX.plainMatrix3D(nabla2[1 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla^2T^i_j = {LaTeX.plainMatrix3D(nabla2[2 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla^3T^i_j = {LaTeX.plainMatrix3D(nabla2[3 - 1])}")}', end=f'{lb}')

            case 4: # nabla^m T^j_i
                doc.write(f'{eq(f"{bs}nabla^1T^j_i = {LaTeX.plainMatrix3D(nabla2[1 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla^2T^j_i = {LaTeX.plainMatrix3D(nabla2[2 - 1])}")}', end=f'{lb}')
                doc.write(f'{eq(f"{bs}nabla^3T^j_i = {LaTeX.plainMatrix3D(nabla2[3 - 1])}")}', end=f'{lb}')

        return [nablaSolutions2, nabla2]


if __name__ == '__main__':
    # packages
    packages = { # from : what
        'inputenc'  : ['utf8'],  
        'babel'     : ['russian'],
        'amsmath'   : [],
        'esdiff'    : ['thinc'],
        'mathtools' : [],
        'fontenc'   : ['T1'],
        'geometry'  : ['top=26mm', 'bottom=20mm']
        # 'geometry'  : ['showframe'],
        # 'layout'    : [],
    }

    # create file for output and LaTeX instance
    with open("latex.tex", "w", encoding="utf-8") as f, LaTeX(f, packages=packages, 
                                                                 title='Ковариантные производные тензорного поля, заданного в цилиндрических координатах.', 
                                                                 author='Очкин Никита Валерьевич ФН11-42Б') as doc:
        
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

        # Tⁱʲ = [[0,   X+Y, 0],
        #        [2*Z, 0,   0],
        #        [0,   0,   0]]

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

                resStr += f' = Q^s_{i} \\cdot Q^p_{j} \\cdot \\delta_{{sp}} = ' 
                
                for s in range(1, 3 + 1):
                    for p in range(1, 3 + 1):

                        mid_res = Qⁱₖ[s - 1][i - 1] * \
                                  Qⁱₖ[p - 1][j - 1] * \
                                    δ[s - 1][p - 1]
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

        QⁱᵐSolutions = [0 for _ in range(3 * 3)]
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
        doc.write(f'По условию')

        TⁱʲSupport = [[0 for _ in range(3)] for _ in range(3)]

        for i in range(1, len(Tⁱʲ) + 1):
            for j in range(1, len(Tⁱʲ[i - 1]) + 1):
                TⁱʲSupport[i - 1][j - 1] = (f'T^{{{i}{j}}} = {LaTeX.cleanUp(Tⁱʲ[i - 1][j - 1])}')

        doc.write(f'{LaTeX.alignat(TⁱʲSupport, quad = True)}')

        # Tᵢⱼ
        T1Solutions = [0 for _ in range(3 * 3)]
        T1          = [[0 for _ in range(3)] for _ in range(3)]

        doc.write(f'Вычислим ковариантные компоненты по формуле:')
        doc.write(f'{eq("T_{ij} = T^{kl}g_{ki}g_{lj};")}')

        count = 0

        for i in range(1, 3 + 1):
            for j in range(1, 3 + 1):

                res = 0
                resStr = f'T_{{{i}{j}}} = '

                resStr += f'T^{{kl}}'
                resStr += rf(f'g_{{k{i}}}g_{{l{j}}} = ')

                for k in range(1, 3 + 1):
                    for l in range(1, 3 + 1):

                        res +=  Tⁱʲ[k - 1][l - 1] * \
                                  g[k - 1][i - 1] * \
                                  g[l - 1][j - 1]

                res     = sp.simplify(res)
                resStr += LaTeX.cleanUp(res)

                T1[i - 1][j - 1]   = res
                T1Solutions[count] = resStr
                count += 1
        
        doc.write(f'{LaTeX.alignat(T1Solutions, quad = True)}')

        doc.write(f'Запишем полученную матрицу:')
        doc.write(f'{eq(f"T_{{ij}}={LaTeX.plainMatrix3D(T1)}")}', end=f'{lb}')

        # Tⁱⱼ
        T2Solutions = [0 for _ in range(3 * 3)]
        T2          = [[0 for _ in range(3)] for _ in range(3)]

        doc.write(f'{eq("T^i_j = T^{ik}g_{kj};")}')

        count = 0

        for i in range(1, 3 + 1):
            for j in range(1, 3 + 1):

                res = 0
                resStr = f'T^{i}_{j} = '

                resStr += f'T^{{{i}{k}}}'
                resStr += rf(f'g_{{k{j}}} = ')

                for k in range(1, 3 + 1):

                    res += Tⁱʲ[i - 1][k - 1] * g[k - 1][j - 1]

                res     = sp.simplify(res)
                resStr += LaTeX.cleanUp(res)

                T2[i - 1][j - 1]   = res
                T2Solutions[count] = resStr
                count += 1
        
        doc.write(f'{LaTeX.alignat(T2Solutions, quad = True)}')

        doc.write(f'Запишем полученную матрицу:')
        doc.write(f'{eq(f"T^i_j={LaTeX.plainMatrix3D(T2)}")}', end=f'{lb}')

        # Tʲᵢ
        T3Solutions = [0 for _ in range(3 * 3)]
        T3          = [[0 for _ in range(3)] for _ in range(3)]

        doc.write(f'{eq("T^j_i = T^{kj}g_{ki};")}')

        count = 0

        for i in range(1, 3 + 1):
            for j in range(1, 3 + 1):

                res = 0
                resStr = f'T^{j}_{i} = '

                resStr += f'T^{{{k}{j}}}'
                resStr += rf(f'g_{{k{i}}} = ')


                for k in range(1, 3 + 1):

                    res += Tⁱʲ[k - 1][j - 1] * g[k - 1][i - 1]

                res     = sp.simplify(res)
                resStr += LaTeX.cleanUp(res)

                T3[i - 1][j - 1]   = res
                T3Solutions[count] = resStr
                count += 1
        
        doc.write(f'{LaTeX.alignat(T3Solutions, quad = True)}')

        doc.write(f'Запишем полученную матрицу:')
        doc.write(f'{eq(f"T^j_i={LaTeX.plainMatrix3D(T3)}")}', end=f'{lb}')

        doc.write(f'Найдем символы Кристоффеля по формуле:')
        doc.write(eq(f"\\Gamma^m_{{ij}} = \\frac{{1}}{{2}}g^{{km}}({bs}frac{{{bs}partial g_{{kj}}}}{{{bs}partial X^i}} + {bs}frac{{{bs}partial g_{{ik}}}}{{{bs}partial X^j}} - {bs}frac{{{bs}partial g_{{ij}}}}{{{bs}partial X^k}}):"))

        # контейнер для символов Кристоффеля 2-го рода
        ГSolutions = [[0 for _ in range(3 * 3)] for _ in range(3)] 
        Г          = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)] # mij

        for m in range(1, 3 + 1):
            count = 0
            for i in range(1, 3 + 1):
                for j in range(1, 3 + 1):
                    
                    res = 0
                    resStr = f'\\Gamma^{m}_{{{i}{j}}} = \\frac{{1}}{{2}}g^{{k{m}}}({bs}frac{{{bs}partial g_{{k{j}}}}}{{{bs}partial X^{i}}} + {bs}frac{{{bs}partial g_{{{i}k}}}}{{{bs}partial X^{j}}} - {bs}frac{{{bs}partial g_{{{i}{j}}}}}{{{bs}partial X^k}}) = '

                    for k in range(1, 3 + 1):
                        res += 1/2 * gInverse[k - 1][m - 1] * (sp.diff(g[k - 1][j - 1], f'X{extra.get_super(i)}') + \
                                                               sp.diff(g[i - 1][k - 1], f'X{extra.get_super(j)}') - \
                                                               sp.diff(g[i - 1][j - 1], f'X{extra.get_super(k)}'))

                    res     = sp.simplify(res)
                    resStr += LaTeX.cleanUp(res)
                    
                    Г[m - 1][i - 1][j - 1]   = res
                    ГSolutions[m - 1][count] = resStr
                    count += 1

        m = LaTeX.m 

        doc.write(f'При m = 1:')
        doc.write(f'{LaTeX.alignat(ГSolutions[1 - 1], quad = True)}')

        doc.write(f'При m = 2:')
        doc.write(f'{LaTeX.alignat(ГSolutions[2 - 1], quad = True)}')

        doc.write(f'При m = 3:')
        doc.write(f'{LaTeX.alignat(ГSolutions[3 - 1], quad = True)}')

        doc.write(f'Запишем результат:')
        doc.write(f'{eq(f"{bs}Gamma^1_{{ij}} = {LaTeX.plainMatrix3D(Г[1 - 1])}")}', end=f'{lb}')
        doc.write(f'{eq(f"{bs}Gamma^2_{{ij}} = {LaTeX.plainMatrix3D(Г[2 - 1])}")}', end=f'{lb}')
        doc.write(f'{eq(f"{bs}Gamma^3_{{ij}} = {LaTeX.plainMatrix3D(Г[3 - 1])}")}', end=f'{lb}')

        # контейнер для ковариантной производной контравариантных компонент
        nablaSolutions1 = [[0 for _ in range(3 * 3)] for _ in range(3)] 
        nabla1          = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)] # nabla_k T^{ij}

        nablaSolutions1, nabla1 = covariant(Tⁱʲ, Г, 1)

        # контейнер для ковариантной производной ковариантных компонент
        nablaSolutions2 = [[0 for _ in range(3 * 3)] for _ in range(3)] 
        nabla2          = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)] # nabla_k T_{ij}

        nablaSolutions2, nabla2 = covariant(T1, Г, 2)

        # контейнер для ковариантной производной смешанных компонент
        nablaSolutions3 = [[0 for _ in range(3 * 3)] for _ in range(3)] 
        nabla3          = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)] # nabla_k T^i_j

        nablaSolutions3, nabla3 = covariant(T2, Г, 3)

        # контейнер для ковариантной производной смешанных компонент
        nablaSolutions4 = [[0 for _ in range(3 * 3)] for _ in range(3)] 
        nabla4          = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)] # nabla_k T^j_i

        nablaSolutions4, nabla4 = covariant(T3, Г, 4)

        doc.write(f'Вычислим контравариантные производные.')

        # контейнер для контравариантной производной от контравариантных компонент
        nablaSolutions5 = [[0 for _ in range(3 * 3)] for _ in range(3)] 
        nabla5          = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)] # nabla^m T^{ij}

        nablaSolutions5, nabla5 = contravariant(gInverse, nabla1, 1)

        # контейнер для контравариантной производной от ковариантных компонент
        nablaSolutions6 = [[0 for _ in range(3 * 3)] for _ in range(3)] 
        nabla6          = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)] # nabla^m T_{ij}

        nablaSolutions6, nabla6 = contravariant(gInverse, nabla2, 2)

        # контейнер для контравариантной производной от смешанных компонент
        nablaSolutions7 = [[0 for _ in range(3 * 3)] for _ in range(3)] 
        nabla7          = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)] # nabla^m T^i_j

        nablaSolutions7, nabla7 = contravariant(gInverse, nabla3, 3)

        # контейнер для контравариантной производной от смешанных компонент
        nablaSolutions8 = [[0 for _ in range(3 * 3)] for _ in range(3)] 
        nabla8          = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)] # nabla^m T^j_i

        nablaSolutions8, nabla8 = contravariant(gInverse, nabla4, 4)
