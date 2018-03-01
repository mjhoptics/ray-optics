import csv
import os


class MapTLA:
    _d = {}

    def __init__(self):
        TLA, CmdFct, IndxQuals, DataType, Quals = range(5)
        if len(self._d) == 0:
            print('initialize dictionary')
            """
            cwd = os.getcwd()
            print(cwd)
            os.chdir(cwd+"/codev")
            print(os.getcwd())
            """
            with open('codev/tla_mapping.csv', 'rU') as f:
                reader = csv.reader(f)
                for row in reader:
                    if row[TLA] is not '':
                        if row[Quals] is not '':
                            row[Quals] = row[Quals].split(',')
                        self._d[row[TLA]] = row[CmdFct:]

    def find(self, tla):
        try:
            return self._d[tla]
        except KeyError:
            return None
