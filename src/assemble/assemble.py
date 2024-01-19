import sys
import glob
import re

to_stdout = True

def iopt_generator(iopt):
    for x in iopt.replace(" ","").split(','):
        if ":" in x:
            a, b = x.split(':')
            for i in range(int(a), int(b)+1):
                yield str(i)
        else:
            yield x

class Csv_info:

    def __init__(self):
        self.csvs = []

    def __call__(self, iopt, lines):
        data = dict(
        angmom = None,
        type_ = None,
        normalized = None,
        angtype = None,
        multiplicity = None,
        npar = None)

        for line in lines:
            for k in data:
                m = re.match(r'^\s*-\s+'+ k.replace("_", "") +'\s*:?=\s*(.+)\s*$', line)
                if m:
                    data[k] = f"{':' if ':=' in line else ''}{m.group(1)}"

        if data['angmom'] is None:
            return
        if data['npar'] is None:
            return
        if data['multiplicity'] is None:
            return

        for iopt_ in iopt_generator(iopt):
            data_copy = data.copy()
            data_copy['iopt'] = iopt_
            for k in data_copy:
                if data_copy[k].startswith(':'):
                    # execute python code, setting local variables
                    data_copy[k] = eval(data_copy[k][1:], dict(iopt=int(iopt_)))

            self.csvs.append(data_copy)

csv_info = Csv_info()

# Load header
with open('makefun_header.f90', 'r') as f:
    header = f.read()

# Load footer
with open('makefun_footer.f90', 'r') as f:
    footer = f.read()

# Load orbitals
# search for all file that fits name pattern orb_*.f90
orb_files = glob.glob('orb_*.f90')
# read files
orbitals = {f.replace("orb_", "").replace(".f90",""): open(f, 'r').read() for f in orb_files}

# Read comment lines in orbitals
for k, v in orbitals.items():
    lines = v.split('\n')
    comment_lines = []
    for i, line in enumerate(lines):
        if line.strip().startswith('!'):
            comment_lines.append(line.strip()[1:])
        else:
            break
    csv_info(k, comment_lines)

# Assemble makefun_out.f90
with open('makefun_out.f90', 'w') as f:
    if to_stdout:
        f = sys.stdout
    f.write(header)
    f.write('select case (iopt)\n')
    for k, v in orbitals.items():
        f.write(f'case ({k})\n')
        f.write(v)
    f.write(footer)

# Assemble makefun_out.csv
header = "iopt,angmom,type_,normalized,angtype,multiplicity,npar"
with open('makefun_out.csv', 'w') as f:
    f.write(f'{header.replace("-","")}\n')
    for csv in csv_info.csvs:
        f.write(','.join([str(csv[x]) for x in header.split(",")]) + '\n')

