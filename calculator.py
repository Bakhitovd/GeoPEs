import streamlit as st
from read_config import read_config
import datetime
from config import default_text
from io import TextIOWrapper

def main():
    st.title('Бахитовка')
    page = st.sidebar.selectbox("Choose a page", ["Настройка параметров",  "Расчет"])
    st.sidebar.title('Параметры')

    if page == "Настройка параметров":
        st.title("Настройка параметров")
        file_config = st.file_uploader('Загрузка файла конфигураций', type='txt', key='file')
        if file_config is not None:
            stringio = TextIOWrapper(file_config)
            current_p = set_config(stringio)
            to_temp_file('temp', current_p)
        else:
            conf = open('temp' + '.txt', 'r')
            current_p = set_config(conf)

        if st.button('Сохранить настройки в файл'):
            to_file(current_p)
    elif page == "Расчет":
        st.header("РАСЧЕТЫ ")
        st.write("Тут будут расчетики ")
        current_p, sh = read_temp_file('temp')
        st.write(current_p)


def set_config(conf):

    p, sh = read_config(conf)

    field = st.sidebar.text_input('Месторождение', 'Field_1', key="field")
    qz = st.sidebar.text_input('ЗАПАСЫ, млн. ст. м3', p['qz'], key="name")
    # PVT
    st.sidebar.write('PVT')
    mdsn = st.sidebar.text_input('Начальн. молярная доля газа сепарации в пластовом газе', p['mdsn'])
    rkon = st.sidebar.text_input('Плотность конденсата, кг/м3', p['rkon'])
    kusk = st.sidebar.text_input('Коэффициент усадки конденсата, д.ед.', p['kusk'])
    mmas = st.sidebar.text_input('Молекулярная масса смеси, г/моль', p['mmas'])
    pot = st.sidebar.text_input('Начальное содержание конденсата, гр/л', p['pot'])
    ppln = st.sidebar.text_input('Давление пластовое начальное, ата', p['ppln'])
    tpl = st.sidebar.text_input('Температура пластовая начальная,К', p['tpl'])
    row = st.sidebar.text_input('Относительная плотность по воздуху ', p['row'])
    ppkr = st.sidebar.text_input('Псевдокритическое давление, ата', p['ppkr'])
    tpkr = st.sidebar.text_input('Псевдокритическая температура, К', p['tpkr'])
    zkr = st.sidebar.text_input('zкр', p['zkr'])
    st.sidebar.write('Кривая конденсато-содержания')
    kp3 = st.sidebar.text_input('k3', list(p['kp'])[0])
    kp2 = st.sidebar.text_input('k2', list(p['kp'])[1])
    kp1 = st.sidebar.text_input('k1', list(p['kp'])[2])
    kp0 = st.sidebar.text_input('k0', list(p['kp'])[3])

    # Обустройство
    st.sidebar.write('Обустройство')
    pbd = st.sidebar.text_input('Минимально допустимое входное Рукпг, ата', p['pbd'])
    lsk = st.sidebar.text_input('Глубина скважин, м', p['lsk'])
    dnkt = st.sidebar.text_input('Диаметр НКТ, мм', p['dnkt'])
    nkust = st.sidebar.text_input('Число скважин в кусте, шт', p['nkust'])
    lsh = st.sidebar.text_input('Длина шлейфа, м', p['lsh'])
    diamsh = st.sidebar.text_input('Диаметр шлейфа, мм', p['diamsh'])
    lkoll = st.sidebar.text_input('Длина коллектора, м', p['lkoll'])
    diamkoll = st.sidebar.text_input('Диаметр коллектора, мм', p['diamkoll'])
    keks = st.sidebar.text_input('Коэффициент эксплуатации', p['keks'])
    ksep = st.sidebar.text_input('Коэффициент сепарации', p['ksep'])
    pesgmin = st.sidebar.text_input('Давление входа в ЕСГ, ата', p['pesgmin'])
    dpukpg = st.sidebar.text_input('Потери давление на УКПГ, ата', p['dpukpg'])

    # Скважина
    st.sidebar.write('Скважина')
    kint = st.sidebar.text_input('Коэффициент интенсификации скважин', p['kint'])
    tust = st.sidebar.text_input('Температура на устье скважин, K', p['tust'])
    deprd = st.sidebar.text_input('Максимально допустимая депрессия, атм', p['deprd'])
    ka = st.sidebar.text_input('коэффициент продуктивности А', p['ka'])
    kb = st.sidebar.text_input('коэффициент продуктивности B', p['kb'])
    qn = st.sidebar.text_input('Начальный средний дебит скважин, тыс.м3/сут', p['qn'])

    # Расписание
    st.sidebar.write('Расписание')
    start_year = st.sidebar.text_input('Начало расчета, год', int(p['start_year']))
    end_year = st.sidebar.text_input('Конец расчета, год', int(p['end_year']))
    st.sidebar.write('Расписание')
    sh_input = st.sidebar.text_area('Ввод скважин, полка добычи млн.м3/год', default_text)


    current_p = {'mdsn': mdsn, 'rkon': rkon, 'kusk': kusk, 'mmas': mmas, 'pot': pot, 'ppln': ppln, 'tpl': tpl, 'row': row, 'ppkr': ppkr, 'tpkr': tpkr, 'zkr': zkr,
                       'kp': kp3 + ' ' + kp2 + ' ' + kp1 + ' ' + kp0,
                       'kint' : kint, 'tust' : tust, 'deprd' : deprd, 'ka' : ka, 'kb' : kb, 'qn' : qn, 'qz' : qz,
                       'pbd' : pbd, 'lsk' : lsk, 'dnkt' : dnkt, 'nkust' : nkust, 'lsh' : lsh, 'diamsh' : diamsh, 'lkoll' : lkoll, 'diamkoll' : diamkoll, 'keks' : keks, 'ksep' : ksep, 'pesgmin' : pesgmin,
                       'dpukpg' : dpukpg, 'start_year' : start_year, 'end_year' : end_year, 'sh_input' : sh_input, 'field': field }
    return current_p


def to_file(paramet):

    today = datetime.datetime.today().strftime("%Y-%m-%d-%H.%M.%S")
    new_conf = open(paramet['field'] + '_config_' + today + '.txt', 'w')
    for key in paramet:
        if key == 'field' or key == 'sh_input':
            continue
        else:
            new_conf.write(key + ' ' + paramet[key] + '\n')
    new_conf.write(paramet['sh_input'])
    new_conf.close()


def read_temp_file(filename):
    conf = open(filename + '.txt', 'r')
    p, sh = read_config(conf)
    return p, sh


def to_temp_file(filename, paramet):

    new_conf = open(filename + '.txt', 'w')
    for key in paramet:
        if key == 'field' or key == 'sh_input':
            continue
        else:
            new_conf.write(key + ' ' + paramet[key] + '\n')
    new_conf.write(paramet['sh_input'])
    new_conf.close()

if __name__ == "__main__":
    main()