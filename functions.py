import math


def C5(kp, ppl):
    """
    Расчитываем конденсато-содержание по палетке
    """
    return kp[0] * ppl ** 3 + kp[1] * ppl ** 2 + kp[2] * ppl + kp[3]


def velosity(param, qn, pz):
    """
    Расчет плотности газа, скорости движения газа на забое и минимально допустимой скорости для выноса жидкой фазы
    """
    rg = rogas(param['ppkr'], param['tpkr'], param['row'], pz, param['tpl'])
    vz = qn * 1.205 * param['row'] / rg / 86.4 / (3.14 * (param['dnkt'] / 1000) ** 2 / 4)
    vmin = 5.73 * (45 - 0.0455 * pz) ** 0.25 * pz ** - 0.5
    return vz, vmin


def kdt(p1, p2, t1):
    # p1 - давление в 1 точке, ата
    # p2 - давление во 2 точке, ата
    # t1 - температурав 1 точке, oC
    # dd1,dd2,dd3 и dd4 - коэф. в уравнении изотальпы, проходящей через  p1 и t1
    # Результат - возвращает температуру во 2 точке, оС
    dd1 = 0.0000309375
    dd2 = 0.023106
    dd3 = - 0.01997775
    dd4 = 5.54827
    fqi = ((math.log(p1) / 2.3026) - dd4 - dd2 * t1) / (dd1 * t1 + dd3)
    t2 = ((math.log(p2) / 2.3026) - dd4 - dd3 * fqi) / (dd1 * fqi + dd2)
    # di=(t1-t2)/(p1-p2)
    fn_return_value = t2
    return fn_return_value


def koefz(tpk, ppk, t, p):
    # tpk - температура критическая, К
    # ppk - давление критическое, ата
    # t - текущая температура, К
    # p - текущее давление, ата
    # Результат - возвращает значение коэффичиента сверхсжимаемости при p и t
    tp = t / tpk
    pp = p / ppk
    fn_return_value = (0.4 * math.log(tp) / 2.303 + 0.73) ** pp + pp / 10
    return fn_return_value


def my(ppkr, tpkr, mol, rg, tsr, psr, ro, zkr):
    # ppkr- критическое давление,  ата
    # tpkr- критическая температура, К
    # mol - молекулярная масса смеси  (для метана ­  16.24)
    # rg  - текущая плотность газа при рsr и tsr, кг/м3
    # tsr - текущая температура, К
    # psr - текущее давление,  ата
    # ro - относительная плотность газа по воздуху
    # zkr - критический коэффициент сверхсжимаемости
    # Результат - возвращает значение вязкости при psr и tsr, Па*с
    e = tpkr ** (1 / 6) / (mol ** 0.5 * ppkr ** (2 / 3))
    # mo = 0.0101 * (tsr - 273.2) ^ (1 / 8) - 0.00576 * ro ^ 0.5
    mo = 0.0101 * tsr ** (1 / 8) - 0.00576 * ro ** 0.5
    # rpr = rg / 162
    z = koefz(tpkr, ppkr, tsr, psr)
    zpr = z / zkr
    ppr = psr / ppkr
    tpr = tsr / tpkr
    rpr = ppr / tpr / zpr
    if psr < 50:
        m = ((
                     0.1023 + 0.023364 * rpr + 0.058533 * rpr ** 2 - 0.040758 * rpr ** 3 + 0.0093324 * rpr ** 4) ** 4 - 0.0001) / e + mo
        m = m / 1000
    else:
        m = (0.000108 * (math.exp(1.439 * rpr) - 1 / math.exp(1.111 * rpr ** 1.858))) / e + mo
        m = m / 1000
    fn_return_value = m
    return fn_return_value


def p3(ppkr, tpkr, ro, mol, pp, tp, py, qgas, tr, v, ty, d, l):
    # ppkr - критическое давление, ата
    # tpkr - критическая температура, К
    # ro   - относ. плотность газа по воздуху
    # mol  - молекулярная масса смеси
    # pp   - пластовое давление, ата
    # tp   - пластовая температура, К
    # py   - устьевое давление, ата
    # qgas - дебит газа, тыс.м3/сут
    # tr   - время режима, мин
    # v    - объем вынесенной жидкости, см3
    # ty   - температура на устье, К
    # d    - диаметр трубы, мм
    # l    - длина трубы, м
    # Результат - возвращает забойное давление в скважине,
    #            расчитанное по движ. столбу смеси в трубе, ата
    tc = tp
    qb = v * 0.000001 / tr * 1.44
    mb = qb * 1000
    pit = py()
    rn = ro * 1.205
    fu = 10000
    while fu > 0.01:
        psr = (py() + pit) / 2
        tsr = (tc - ty) / math.log(tc / ty)
        kzsr = koefz(tpkr, ppkr, tsr, psr)
        rg = rn * psr * 293 / (1.033 * kzsr * tsr)
        qr = qgas * rn / rg
        fi = qr / (qr + qb)
        r = fi + (1 - fi) * 1000 / rg
        s = 0.03415 * ro * r * l / (kzsr * tsr)
        dcm = d / 10
        her = 2 * 0.12 / (10 * dcm)
        gg = qgas * rg
        qsm = (gg + mb) / rg
        e = tpkr ** (1 / 6) / (mol ** 0.5 * ppkr ** (2 / 3))
        # mo = 0.0101 * (tsr - 273) ^ (1 / 8) - 0.00576 * ro ^ 0.5
        mo = 0.0101 * tsr ** (1 / 8) - 0.00576 * ro ** 0.5
        rpr = rg / 162
        if psr < 50:
            m = ((
                         0.1023 + 0.023364 * rpr + 0.058533 * rpr ** 2 - 0.040758 * rpr ** 3 + 0.0093324 * rpr ** 4) ** 4 - 0.0001) / e + mo
            m = m / 1000
        else:
            m = (0.000108 * (math.exp(1.439 * rpr) - 1 / math.exp(1.111 * rpr ** 1.858))) / e + mo
            m = m / 1000
        qgrc = qr / 86.4
        dmt = d / 1000
        vgrc = qgrc / (3.14 * dmt ** 2 / 4)
        prey = vgrc * dmt * rg / m
        # lq=(1/(2*log10(7.41/her)))^2
        lq = 1 / (1.8 * math.log(prey) / 2.303 - 1.64) ** 2
        if prey < 100000:
            lq = 0.3164 / (prey) ** 0.25
        p3a = (py() ** 2 * math.exp(2 * s) + 1.377 * lq * kzsr ** 2 * tsr ** 2 * qsm ** 2 / (r * dcm ** 5) * (
                math.exp(2 * s) - 1)) ** 0.5
        fu = abs(pit - p3a)
        pit = p3a
        # tc=fnkdt(pp,p3,tp-273.2)+273.2
        # ? lq,prey,her,qsm,p3
    fn_return_value = p3a
    return fn_return_value


def pbar(ppkr, tpkr, ro, tp, pct, l, tct):
    # ppkr - критическое давление, ата
    # tpkr - критическая температура, К
    # ro   - относительная плотность газа по воздуху
    # tp   - пластовая температура, К
    # pct  - статическое давление на устье, ата
    # l    - длина трубы, м
    # tct  - температура на устье, К
    # Результат - возвращает забойное давление в скважине,
    #            расчитанное по неподвижному столбу газа (барометрическая формула), ата
    pit = pct
    fu = 1000
    tsr = (tp - tct) / math.log(tp / tct)
    while fu > 0.01:
        psr = (pct + pit) / 2
        kzsr = koefz(tpkr, ppkr, tsr, psr)
        s = 0.03415 * ro * l / (kzsr * tsr)
        ppl = pct * math.exp(s)
        fu = abs(ppl - pit)
        pit = ppl
        # PRINT ppl, kzsr, s, l
        # WHILE INKEY$ = "": WEND
    fn_return_value = ppl
    return fn_return_value


def q1(r, d, t, p1, p2, diam, tpk, ppk):
    # r - относительная плотность газа по воздуху
    # d - диаметр диафрагмы, мм
    # t - температура перед диафрагмой, К
    # p1- давление перед диафрагмой, ата
    # p2-давление после диафрагмы, ата
    # diam - диаметр дикта или трубы, мм
    # tpk - критическая температура, К
    # ppk - критическое давление, ата
    # Результат - возвращает замеренный дебит газа, тыс.м3/сут
    zz = koefz(tpk, ppk, t, p1)
    select_variable_0 = p2
    if (select_variable_0 > 0):
        dt = diam / 1000
        dd = d / 1000
        f0 = 3.14 * dd ** 2 / 4
        f1 = 3.14 * dt ** 2 / 4
        m = f0 / f1
        bet = 1 / (1 - m ** 2) ** 0.5
        al = 0.61 * bet
        popr = 1 - (0.3707 + 0.3184 * m ** 2) * (1 - (p2 / p1) ** (1 / 1.35)) ** 0.935
        rn = r * 1.205
        rg = rn * p1 * 293 / (1.033 * t * zz)
        q = al * popr * f0 * Sqr(200000 * (p1 - p2) / rg)
        q = q * 86.4 * rg / rn
    elif (select_variable_0 == 0):
        fn_return_value = (0.183 * d ** 2 * p1) / Sqr(r * zz * t)
        if diam < 90:
            fn_return_value = (0.1802 * d ** 2 * p1) / Sqr(r * zz * t)
    return fn_return_value


def qab1(a, B, c, ppl, pzab):
    # a, b, c   - коэффициенты в формуле притока
    # ppl - пластовое давление, ата
    # pzab- забойное давление, ата
    # Результат - возвращает дебит скважины при заданных ppl и pzab
    fun = 1000000
    qmin = - 100000
    qmax = 100000
    iter = 0
    while abs(fun) > 0.001:
        iter = iter + 1
        qi = (qmax + qmin) / 2
        fun = a * qi + B * qi ** 2 + c - ppl ** 2 + pzab ** 2
        if fun < 0:
            qmin = qi
        else:
            qmax = qi
        if iter > 2000:
            qi = 0
            break
        # PRINT qi, fun
        # WHILE INKEY$ = "": WEND
    # PRINT "------------"
    fn_return_value = qi
    return fn_return_value


def rogas(ppk, tpk, ro, psr, tsr):
    # ppk - критическое давление, ата
    # tpk - критическая температура, К
    # ro  - относ. плотность газа по воздуху
    # psr - текущее давление, ата
    # tsr - текущая температура, К
    # Результат - возвращает плотность газа при psr и tsr
    rn = 1.205 * ro
    zsr = koefz(tpk, ppk, tsr, psr)
    rg = rn * psr * 293 / (1.033 * zsr * tsr)
    fn_return_value = rg
    # print(ppk, tpk, ro, psr, tsr)
    return fn_return_value


def pyz(ppkr, tpkr, ro, mol, tp, pz, qgas, tr, v, ty, d, l, rl):
    tsr = 0

    psr = 0
    # ppkr - критическое давление, ата
    # tpkr - критическая температура, К
    # ro   - относ. плотность газа по воздуху
    # mol  - молекулярная масса смеси
    # tp   - пластовая температура, К
    # pz   - забойное давление, ата
    # qgas - дебит газа, тыс.м3/сут
    # tr   - время режима, мин
    # v    - объем вынесенной жидкости, см3
    # ty   - температура на устье, К
    # d    - диаметр трубы, мм
    # l    - длина трубы, м
    # rl - плотность жидкости, кг/м3
    # Результат - возвращает устьевое давление в скважине,
    #            расчитанное по движ. столбу смеси в трубе, ата
    # VB2PY (UntranslatedCode) On Error GoTo 10
    tc = tp
    qb = v * 0.000001 / tr * 1.44
    mb = qb * rl
    pit = pz
    rn = ro * 1.205
    fu = 10000
    if qgas == 0:
        fn_return_value = 0
        return fn_return_value
    while fu > 0.0001:
        psr = (pz + pit) / 2
        tsr = (tc - ty) / math.log(tc / ty)
        kzsr = koefz(tpkr, ppkr, tsr, psr)
        rg = rn * psr * 293 / (1.033 * kzsr * tsr)
        qr = qgas * rn / rg
        fi = qr / (qr + qb)
        r = fi + (1 - fi) * rl / rg
        s = 0.03415 * ro * r * l / (kzsr * tsr)
        dcm = d / 10
        her = 2 * 0.12 / (10 * dcm)
        gg = qgas * rg
        qsm = (gg + mb) / rg
        e = tpkr ** (1 / 6) / (mol ** 0.5 * ppkr ** (2 / 3))
        # mo = 0.0101 * (tsr - 273) ^ (1 / 8) - 0.00576 * ro ^ 0.5
        mo = 0.0101 * tsr ** (1 / 8) - 0.00576 * ro ** 0.5
        rpr = rg / 162
        if psr < 50:
            m = ((
                         0.1023 + 0.023364 * rpr + 0.058533 * rpr ** 2 - 0.040758 * rpr ** 3 + 0.0093324 * rpr ** 4) ** 4 - 0.0001) / e + mo
            m = m / 1000
        else:
            m = (0.000108 * (math.exp(1.439 * rpr) - 1 / math.exp(1.111 * rpr ** 1.858))) / e + mo
            m = m / 1000
        qgrc = qr / 86.4
        dmt = d / 1000
        vgrc = qgrc / (3.14 * dmt ** 2 / 4)
        prey = vgrc * dmt * rg / m
        # lq=(1/(2*log10(7.41/her)))^2
        lq = 1 / (1.8 * math.log(prey) / 2.303 - 1.64) ** 2
        if prey < 100000:
            lq = 0.3164 / (prey) ** 0.25
        lq = 0.017
        xxx = pz ** 2
        yyy = 1.377 * lq * kzsr ** 2 * tsr ** 2 * qsm ** 2 / (r * dcm ** 5) * (math.exp(2 * s) - 1)
        p3a = ((pz ** 2 + 1.377 * lq * kzsr ** 2 * tsr ** 2 * qsm ** 2 / (r * dcm ** 5) * (
                math.exp(2 * s) - 1)) / math.exp(2 * s)) ** 0.5
        # p3a = (pz ^ 2 - 1.377 * lq * kzsr ^ 2 * tsr ^ 2 * qgas ^ 2 / dcm ^ 5 * (math.exp(2 * s) - 1)) ^ 0.5
        fu = abs(pit - p3a)
        pit = p3a
        # tc=fnkdt(pp,p3,tp-273.2)+273.2
        # ? lq,prey,her,qsm,p3
    fn_return_value = p3a
    return fn_return_value
    fn_return_value = - 1
    return fn_return_value


def pipe(pkr, tkr, roo, p1, t, l, diam, doby):
    ptek = 0.0
    # /pkr = 46
    # tkr = 212
    # roo = 0.6944
    # p1 = Лист1.Cells(2, 2).Value
    # t = Лист1.Cells(2, 3).Value
    # l = Лист1.Cells(2, 5).Value
    # diam = Лист1.Cells(2, 6).Value
    # doby = Лист1.Cells(2, 7).Value
    lq = 0.017
    if p1 == - 1:
        fn_return_value = - 1
        return fn_return_value
    dobysek = doby * 1000 / 24 / 60 / 60
    dobysut = doby * 1000000000 / 365
    deltal = 1
    lras = 0
    dps = 0
    deltap = 0
    ptek = p1
    dm = diam / 1000
    plo = 3.14 * dm ** 2 / 4
    while lras <= l:
        lras = lras + deltal
        ptek = ptek - deltap
        if ptek < 0:
            ptek = - 1
            break

        rog = rogas(pkr, tkr, roo, ptek, t)
        # print(rog)
        qreal = dobysek * 1.205 * roo / rog
        w = qreal / plo
        deltap = 0.5 * lq * rog * deltal * w ** 2 / dm / 100000
    p2 = ptek
    # Лист1.Cells(2, 4).Value = rog
    # Лист1.Cells(2, 8).Value = p2
    fn_return_value = p2
    if p2 < 0:
        fn_return_value = - 1
        return fn_return_value
    return fn_return_value


def pipe1(pkr, tkr, roo, p1, t, l, diam, doby):
    ptek = 0
    # /pkr = 46
    # tkr = 212
    # roo = 0.6944
    # p1 = Лист1.Cells(2, 2).Value
    # t = Лист1.Cells(2, 3).Value
    # l = Лист1.Cells(2, 5).Value
    # diam = Лист1.Cells(2, 6).Value
    # doby = Лист1.Cells(2, 7).Value
    lq = 0.017
    if p1 == - 1:
        fn_return_value = - 1
        return fn_return_value
    dobysek = doby * 1000 / 24 / 60 / 60
    dobysut = doby * 1000000000 / 365
    deltal = 1
    lras = 0
    dps = 0
    deltap = 0
    ptek = p1
    dm = diam / 1000
    plo = 3.14 * dm ** 2 / 4
    while lras <= l:
        lras = lras + deltal
        ptek = ptek + deltap
        if ptek < 0:
            ptek = - 1
            break
        rog = rogas(pkr, tkr, roo, ptek, t)
        qreal = dobysek * 1.205 * roo / rog
        w = qreal / plo
        deltap = 0.5 * lq * rog * deltal * w ** 2 / dm / 100000
    p2 = ptek
    # Лист1.Cells(2, 4).Value = rog
    # Лист1.Cells(2, 8).Value = p2
    fn_return_value = p2
    if p2 < 0:
        fn_return_value = - 1
        return fn_return_value
    return fn_return_value


def zzz(pkr, tkr, p, t):
    import math
    # pkr - псевдокритическое давление, ата
    # tkr - псевдокритическая температура, К
    # p   - текущее давление, ата
    # t   - текущая температура, К
    # Результат - возвращает значение коэфф. сверхсжимаемости
    fn_return_value = (0.4 * math.log(t / tkr) / 2.303 + 0.73) ** (p / pkr) + 0.1 * (p / pkr)
    return fn_return_value


def lll(a):
    fn_return_value = math.log(a)
    return fn_return_value


# ppl(ppln, pplsr[i - 1], zn, qz, qsg[i] - qszak[i], ppkr, tpkr, tpl)

def ppl(ppn, pp, zn, qz, qd, pkr1, tkr1, tpl):
    z = 0.0

    iter1 = 0.0

    iter2 = 0.0
    # ppn - начальное пластовое давление, ата
    # pp - пластовое давление за предыдщий год, ата
    # zn - z при начальном пластовом давлении
    # qz - запасы газа, млн.м3
    # qd - суммарное добытое количество газа, млн.м3
    # pkr1 - критическое давление, ата
    # tkr1 - критическая температура, К
    # tpl - пластовая температура, К
    # Результат - возвращает пластовое давление по мат. балансу, ата
    iter2 = pp
    z = zzz(pkr1, tkr1, pp, tpl)
    # Лист2.Cells(1, 3).Value = pp
    iter1 = ppn * z * (1 - qd / qz) / zn
    j = 0
    while abs(iter1 - iter2) > 0.1:
        j = j + 1
        iter2 = iter1
        z = zzz(pkr1, tkr1, iter1, tpl)
        iter1 = ppn * z * (1 - qd / qz) / zn
    # Лист2.Cells(1, 5).Value = j
    fn_return_value = iter1
    return fn_return_value


def qab(a, B, c, ppl, pzab):
    fun = 0.0

    iter = 0

    qi = 0.0
    # a, b, c   - коэффициенты в формуле притока
    # ppl - пластовое давление, ата
    # pzab- забойное давление, ата
    # Результат - возвращает дебит скважины при заданных ppl и pzab
    fun = 1000000
    qmin = - 100000
    qmax = 100000
    iter = 0
    while abs(fun) > 0.1:
        iter = iter + 1
        qi = (qmax + qmin) / 2
        fun = a * qi + B * qi ** 2 + c - ppl ** 2 + pzab ** 2
        if fun < 0:
            qmin = qi
        else:
            qmax = qi
        if iter > 2000:
            qi = 0
            break
        # PRINT qi, fun
        # WHILE INKEY$ = "": WEND
    # PRINT "------------"
    fn_return_value = qi
    return fn_return_value


def py(ppkr, tpkr, ro, mol, tp, pz, qgas, tr, v, ty, d, l, rl):
    tsr = 0.0

    psr = 0.0
    # ppkr - критическое давление, ата
    # tpkr - критическая температура, К
    # ro   - относ. плотность газа по воздуху
    # mol  - молекулярная масса смеси
    # tp   - пластовая температура, К
    # pz   - забойное давление, ата
    # qgas - дебит газа, тыс.м3/сут
    # tr   - время режима, мин
    # v    - объем вынесенной жидкости, см3
    # ty   - температура на устье, К
    # d    - диаметр трубы, мм
    # l    - длина трубы, м
    # rl - плотность жидкости, кг/м3
    # Результат - возвращает устьевое давление в скважине,
    #            расчитанное по движ. столбу смеси в трубе, ата

    # VB2PY (UntranslatedCode) On Error GoTo 10
    try:
        tc = tp
        qb = v * 0.000001 / tr * 1.44
        mb = qb * rl
        pit = pz
        rn = ro * 1.205
        fu = 10000
        if qgas == 0:
            fn_return_value = 0
            return fn_return_value
        while fu > 0.0001:
            psr = (pz + pit) / 2
            tsr = (tc - ty) / math.log(tc / ty)
            kzsr = koefz(tpkr, ppkr, tsr, psr)
            rg = rn * psr * 293 / (1.033 * kzsr * tsr)
            qr = qgas * rn / rg
            fi = qr / (qr + qb)
            r = fi + (1 - fi) * rl / rg
            s = 0.03415 * ro * r * l / (kzsr * tsr)
            dcm = d / 10
            her = 2 * 0.12 / (10 * dcm)
            gg = qgas * rg
            qsm = (gg + mb) / rg
            e = tpkr ** (1 / 6) / (mol ** 0.5 * ppkr ** (2 / 3))
            # mo = 0.0101 * (tsr - 273) ^ (1 / 8) - 0.00576 * ro ^ 0.5
            mo = 0.0101 * tsr ** (1 / 8) - 0.00576 * ro ** 0.5
            rpr = rg / 162
            # print('pz ',pz,'pit ' ,pit)
            if psr < 50:
                m = ((
                             0.1023 + 0.023364 * rpr + 0.058533 * rpr ** 2 - 0.040758 * rpr ** 3 + 0.0093324 * rpr ** 4) ** 4 - 0.0001) / e + mo
                m = m / 1000
            else:
                m = (0.000108 * (math.exp(1.439 * rpr) - 1 / math.exp(1.111 * rpr ** 1.858))) / e + mo
                m = m / 1000
            qgrc = qr / 86.4
            dmt = d / 1000
            vgrc = qgrc / (3.14 * dmt ** 2 / 4)
            prey = vgrc * dmt * rg / m
            # lq=(1/(2*log10(7.41/her)))^2
            lq = 1 / (1.8 * math.log(prey) / 2.303 - 1.64) ** 2
            if prey < 100000:
                lq = 0.3164 / (prey) ** 0.25
            lq = 0.017
            xxx = pz ** 2
            yyy = 1.377 * lq * kzsr ** 2 * tsr ** 2 * qsm ** 2 / (r * dcm ** 5) * (math.exp(2 * s) - 1)
            p3a = ((pz ** 2 - 1.377 * lq * kzsr ** 2 * tsr ** 2 * qsm ** 2 / (r * dcm ** 5) * (
                    math.exp(2 * s) - 1)) / math.exp(2 * s)) ** 0.5
            #  print('params ','pz',pz,lq,kzsr,tsr,'qsm',qsm,r,dcm, s)
            # print('p3a ', p3a)
            # p3a = (pz ^ 2 - 1.377 * lq * kzsr ^ 2 * tsr ^ 2 * qgas ^ 2 / dcm ^ 5 * (math.exp(2 * s) - 1)) ^ 0.5
            fu = abs(pit - p3a)
            pit = p3a
            # tc=fnkdt(pp,p3,tp-273.2)+273.2
            # ? lq,prey,her,qsm,p3
        fn_return_value = p3a
        return fn_return_value
    except:
        fn_return_value = - 1
        return fn_return_value
