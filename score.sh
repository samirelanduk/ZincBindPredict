rm data/csv/*/*.csv

python data/generate_sequence_data.py

python -c "print(); import os; \
    [print(f.split('.')[0].rjust(4), len(open('data/csv/sequence/' + f).readlines()))\
    for f in reversed(sorted(\
        os.listdir('data/csv/sequence'),\
        key=lambda x: len(open('data/csv/sequence/' + x).readlines())\
    )) if f[0] != '.']"

python predict/build_models.py --categories=sequence > /dev/null 2>&1

python -c "print(); import os; import joblib; import warnings; warnings.warn = lambda *a, **k: None;\
    [print(f.split('.')[0].rjust(4), round(joblib.load('predict/models/sequence/' + f.split('.')[0] + '-RF.joblib').matthew_, 2))\
    for f in reversed(sorted(\
        os.listdir('data/csv/sequence'),\
        key=lambda x: len(open('data/csv/sequence/' + x).readlines())\
    )) if f[0] != '.']; print()"

python predict/sequence_tests.py