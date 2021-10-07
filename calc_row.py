import functions as gp
class calc_rows:
    '''
    дебит с учетом коэффициента интенсификации
    '''

    def __init__(self, year, param, ppl, nsk, qg_nak, qk_nak, qz_nak, q_nak_sep):
        self.year = int(year)
        self.param = param
        self.ppl = ppl
        self.nsk = int(nsk)
        self.qg_nak = qg_nak
        self.qk_nak = qk_nak
        self.qz_nak = qz_nak
        self.q_nak_sep = q_nak_sep

    def const_q(self, qn, max_depr):
        '''
        дебит с учетом коэффициента интенсификации
        '''
        self.qn = qn
        self.max_depr = max_depr
        msr = my = gp.my(self.param['ppkr'], self.param['tpkr'], self.param['mmas'], 0, self.param['tpl'], self.ppl,
                         self.param['row'], self.param['zkr']) * 1000
        zsr = gp.zzz(self.param['ppkr'], self.param['tpkr'], self.ppl, self.param['tpl'])
        depr = self.ppl - (self.ppl ** 2 - self.param['ka'] * (self.qn / self.param['kint']) - self.param['kb'] * (
                    self.qn / self.param['kint']) ** 2) ** 0.5
        if (depr.imag != 0):
            return 1
        if depr > self.max_depr:
            return 2
        qg_year = self.nsk * self.qn * 365 * self.param['keks'] / 1000
        qg_nak = self.qg_nak + qg_year
        zn = gp.zzz(self.param['ppkr'], self.param['tpkr'], self.param['ppln'], self.param['tpl'])
        ppl_next = gp.ppl(self.param['ppln'], self.ppl, zn, self.param['qz'], qg_nak, self.param['ppkr'],
                          self.param['tpkr'], self.param['tpl'])
        pz = self.ppl - depr
        pb = gp.py(self.param['ppkr'], self.param['tpkr'], self.param['row'], self.param['mmas'], self.param['tpl'],
                   pz, self.qn, 1440, 0, self.param['tust'], self.param['dnkt'], self.param['lsk'],
                   self.param['rkon'])  # vl добавить потом
        puk = gp.pipe(self.param['ppkr'], self.param['tpkr'], self.param['row'], pb, self.param['tust'],
                      self.param['lsh'],
                      self.param['diamsh'], self.param['nkust'] * self.qn)
        puk = puk - (pb - puk) * 0.15
        pot = gp.C5(self.param['kp'], self.ppl)
        q_year_sep = qg_year * self.param['mdsn']
        q_nak_sep = self.q_nak_sep + q_year_sep
        qk_year = q_year_sep * pot * self.param['ksep'] / 1000
        qk_nak = self.qk_nak + qk_year
        kik = qk_nak / (self.param['pot'] * self.param['qz'] / 1000)
        kig = q_nak_sep / self.param['qz']
        vz, vmin = gp.velosity(self.param, self.qn, pz)
        return {'year' : self.year,
                'nsk': self.nsk,
                'qn': self.qn,
                'depr': depr,
                'ppl': ppl_next,
                'qg_year': q_year_sep,
                'qg_nak_geo': qg_nak,
                'q_nak_sep': q_nak_sep,
                'pot': pot,
                'qk_year': qk_year,
                'qk_nak': qk_nak,
                'kik': kik,
                'pz': pz,
                'pb': pb,
                'puk': puk,
                'vz': vz,
                'vmin': vmin,
                'kig': kig,
                'kik': kik

                }

    def const_depr(self, depr):
        self.depr = depr
        msr = my = gp.my(self.param['ppkr'], self.param['tpkr'], self.param['mmas'], 0, self.param['tpl'], self.ppl,
                         self.param['row'], self.param['zkr']) * 1000
        zsr = gp.zzz(self.param['ppkr'], self.param['tpkr'], self.ppl, self.param['tpl'])
        pz = self.ppl - self.depr
        qn = self.param['kint'] * gp.qab(self.param['ka'], self.param['kb'], 0, self.ppl, pz)
        qg_year = self.nsk * qn * 365 * self.param['keks'] / 1000
        qg_nak = self.qg_nak + qg_year
        zn = gp.zzz(self.param['ppkr'], self.param['tpkr'], self.param['ppln'], self.param['tpl'])
        ppl_next = gp.ppl(self.param['ppln'], self.ppl, zn, self.param['qz'], qg_nak, self.param['ppkr'],
                          self.param['tpkr'], self.param['tpl'])
        pb = gp.py(self.param['ppkr'], self.param['tpkr'], self.param['row'], self.param['mmas'], self.param['tpl'],
                   pz, qn, 1440, 0, self.param['tust'], self.param['dnkt'], self.param['lsk'],
                   self.param['rkon'])  # vl добавить потом
        puk = gp.pipe(self.param['ppkr'], self.param['tpkr'], self.param['row'], pb, self.param['tust'],
                      self.param['lsh'],
                      self.param['diamsh'], self.param['nkust'] * qn)
        puk = puk - (pb - puk) * 0.15
        if puk < self.param['pbd']:
            return 1
        pot = gp.C5(self.param['kp'], self.ppl)
        q_year_sep = qg_year * self.param['mdsn']
        q_nak_sep = self.q_nak_sep + q_year_sep
        qk_year = q_year_sep * pot * self.param['ksep'] / 1000
        qk_nak = self.qk_nak + qk_year
        kik = qk_nak / (self.param['pot'] * self.param['qz'] / 1000)
        kig = q_nak_sep / self.param['qz']
        vz, vmin = gp.velosity(self.param, qn, pz)
        return {'year' : self.year,
                'nsk': self.nsk,
                'qn': qn,
                'depr': self.depr,
                'ppl': ppl_next,
                'qg_year': q_year_sep,
                'qg_nak_geo': qg_nak,
                'q_nak_sep': q_nak_sep,
                'pot': pot,
                'qk_year': qk_year,
                'qk_nak': qk_nak,
                'kik': kik,
                'pz': pz,
                'pb': pb,
                'puk': puk,
                'vz': vz,
                'vmin': vmin,
                'kig': kig,
                'kik': kik
                }

    def decline_production(self):

        msr = my = gp.my(self.param['ppkr'], self.param['tpkr'], self.param['mmas'], 0, self.param['tpl'], self.ppl,
                         self.param['row'], self.param['zkr']) * 1000
        zsr = gp.zzz(self.param['ppkr'], self.param['tpkr'], self.ppl, self.param['tpl'])

        fun = 10000
        deprmin = 0
        deprmax = self.ppl - 1
        while abs(fun) > 0.01:
            depress = (deprmin + deprmax) / 2
            pz = self.ppl - depress
            qn = self.param['kint'] * gp.qab(self.param['ka'], self.param['kb'], 0, self.ppl, pz)

            pb = gp.py(self.param['ppkr'], self.param['tpkr'], self.param['row'], self.param['mmas'], self.param['tpl'],
                       pz, qn, 1440, 0, self.param['tust'], self.param['dnkt'], self.param['lsk'],
                       self.param['rkon'])  # vl добавить потом
            puk = gp.pipe(self.param['ppkr'], self.param['tpkr'], self.param['row'], pb, self.param['tust'],
                          self.param['lsh'],
                          self.param['diamsh'], self.param['nkust'] * qn)
            puk = puk - (pb - puk) * 0.15
            fun = puk - self.param['pbd']
            if puk == - 1:
                fun = - 1000
            if fun < 0:
                deprmax = depress
            else:
                deprmin = depress
        qg_year = self.nsk * qn * 365 * self.param['keks'] / 1000
        qg_nak = self.qg_nak + qg_year

        zn = gp.zzz(self.param['ppkr'], self.param['tpkr'], self.param['ppln'], self.param['tpl'])
        ppl_next = gp.ppl(self.param['ppln'], self.ppl, zn, self.param['qz'], qg_nak, self.param['ppkr'],
                          self.param['tpkr'], self.param['tpl'])
        pb = gp.py(self.param['ppkr'], self.param['tpkr'], self.param['row'], self.param['mmas'], self.param['tpl'],
                   pz, qn, 1440, 0, self.param['tust'], self.param['dnkt'], self.param['lsk'],
                   self.param['rkon'])  # vl добавить потом
        puk = gp.pipe(self.param['ppkr'], self.param['tpkr'], self.param['row'], pb, self.param['tust'],
                      self.param['lsh'],
                      self.param['diamsh'], self.param['nkust'] * qn)
        puk = puk - (pb - puk) * 0.15
        pot = gp.C5(self.param['kp'], self.ppl)
        q_year_sep = qg_year * self.param['mdsn']
        q_nak_sep = self.q_nak_sep + q_year_sep
        qk_year = q_year_sep * pot * self.param['ksep'] / 1000
        qk_nak = self.qk_nak + qk_year
        kik = qk_nak / (self.param['pot'] * self.param['qz'] / 1000)
        kig = q_nak_sep / self.param['qz']
        vz, vmin = gp.velosity(self.param, qn, pz)
        return {'year' : self.year,
                'nsk': self.nsk,
                'qn': qn,
                'depr': depress,
                'ppl': ppl_next,
                'qg_year': q_year_sep,
                'qg_nak_geo': qg_nak,
                'q_nak_sep': q_nak_sep,
                'pot': pot,
                'qk_year': qk_year,
                'qk_nak': qk_nak,
                'kik': kik,
                'pz': pz,
                'pb': pb,
                'puk': puk,
                'vz': vz,
                'vmin': vmin,
                'kig': kig,
                'kik': kik
                }