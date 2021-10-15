def main_func(con):
    import pandas as pd
    from calc_row import calc_rows
    from read_config import read_config

    p, sh = read_config(con)
    deprd = p['deprd']
    ppl = p['ppln']
    qg_nak, qk_nak, q_nak_sep, qz_nak, flag_1, vmin, vz = 0, 0, 0, 0, 0, 0, 0
    pbd = p['keks']
    df = pd.DataFrame({'year': [], 'nsk': [], 'qn': [], 'depr': [], 'ppl': [], 'qg_year': [], 'qg_nak_geo': [], 'pz': [], 'pb': [], 'puk': []})

    for i, row in sh.iterrows():
        year = int(row['year'])
        nsk = int(row['nsk'])
        Q_max = float(row['qmax']) / 365
        qn = Q_max / p['keks'] / nsk
        row = {}
        c = calc_rows(year, p, ppl, nsk, qg_nak, qk_nak, qz_nak, q_nak_sep)
        if flag_1 == 0:
            row_q = c.const_q(qn, deprd)
            row = row_q
            if row_q == 1 or row_q == 2:
                row_depr = c.const_depr(deprd)
                row = row_depr
                if row_depr != 1:
                    df = df.append(row_depr, ignore_index=True)
                    ppl = row_depr['ppl']
                    qg_nak = row_depr['qg_nak_geo']
                    q_nak_sep = row_depr['q_nak_sep']
                    qk_nak = row_depr['qk_nak']
                    vz = row_depr['vz']
                    vmin = row_depr['vmin']
                else:
                    flag_1 = 1
            else:
                df = df.append(row_q, ignore_index=True)
                ppl = row_q['ppl']
                qg_nak = row_q['qg_nak_geo']
                q_nak_sep = row_q['q_nak_sep']
                qk_nak = row_q['qk_nak']
                vz = row_q['vz']
                vmin = row_q['vmin']
        else:
            row_decl = c.decline_production()
            row = row_decl
            df = df.append(row_decl, ignore_index=True)
            ppl = row_decl['ppl']
            qg_nak = row_decl['qg_nak_geo']
            q_nak_sep = row_decl['q_nak_sep']
            qk_nak = row_decl['qk_nak']
            vz = row_decl['vz']
            vmin = row_decl['vmin']
        if vz < vmin:
            break
    df[['year', 'nsk']] = df[['year', 'nsk']].astype(int)
    df[['qn', 'depr', 'ppl', 'qg_year', 'qg_nak_geo', 'pz', 'pb', 'puk']] = df[['qn', 'depr', 'ppl', 'qg_year', 'qg_nak_geo', 'pz', 'pb', 'puk']].astype(float)
    return df


if __name__ == "__main__":
    conf = open('config.txt', 'r')
    print(main_func(conf))