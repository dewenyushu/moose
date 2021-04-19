import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv

x_coord = []
contact_pressure = []


fig1, ax1 = plt.subplots()
# fig2, ax2 = plt.subplots()

# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# nodface
x_coord = []
contact_pressure = []
with open('antonio_chkfile_contact_post_0030.csv') as csv_file: #293
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            # print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            contact_pressure.append(float(row[0]))
            x_coord.append(float(row[2]))
            # print(f'\t{row[0]} works in the {row[1]} department, and was born in {row[2]}.')
            line_count += 1
    print(f'Processed {line_count} lines.')

contact_pressure_ = [x for _,x in sorted(zip(x_coord,contact_pressure))]
contact_pressure = contact_pressure_
# XXX: .sort()
ax1.plot(x_coord,contact_pressure, 'b', linewidth=2, label='Node-on-segment')

x_coord = []
contact_pressure = []
with open('antonio.standard_chkfile_cont_press_0030.csv') as csv_file: #293
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            # print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            contact_pressure.append(float(row[0]))
            x_coord.append(float(row[2]))
            # print(f'\t{row[0]} works in the {row[1]} department, and was born in {row[2]}.')
            line_count += 1
    print(f'Processed {line_count} lines.')

contact_pressure_ = [x for _,x in sorted(zip(x_coord, contact_pressure))]
contact_pressure = contact_pressure_
x_coord.sort()
ax1.plot(x_coord,contact_pressure, 'r--', linewidth=2, label='Standard mortar')

x_coord = []
contact_pressure = []
with open('antonio.standard.second_chkfile_cont_press_0030.csv') as csv_file: #293
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            # print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            contact_pressure.append(float(row[0]))
            x_coord.append(float(row[2]))
            # print(f'\t{row[0]} works in the {row[1]} department, and was born in {row[2]}.')
            line_count += 1
    print(f'Processed {line_count} lines.')

contact_pressure_ = [x for _,x in sorted(zip(x_coord, contact_pressure))]
contact_pressure = contact_pressure_
x_coord.sort()
ax1.plot(x_coord,contact_pressure, 'y--', linewidth=2, label='Standard mortar mixed')

x_coord = []
contact_pressure = []
with open('antonio.standard.second.second_chkfile_cont_press_0030.csv') as csv_file: #293
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            # print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            contact_pressure.append(float(row[0]))
            x_coord.append(float(row[2]))
            # print(f'\t{row[0]} works in the {row[1]} department, and was born in {row[2]}.')
            line_count += 1
    print(f'Processed {line_count} lines.')

contact_pressure_ = [x for _,x in sorted(zip(x_coord, contact_pressure))]
contact_pressure = contact_pressure_
x_coord.sort()
ax1.plot(x_coord,contact_pressure, 'k--', linewidth=2, label='Standard mortar 2nd-2nd')

x_coord = []
contact_pressure = []
with open('antonio.dual_chkfile_cont_press_0030.csv') as csv_file: #293
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            # print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            contact_pressure.append(float(row[0]))
            x_coord.append(float(row[2]))
            # print(f'\t{row[0]} works in the {row[1]} department, and was born in {row[2]}.')
            line_count += 1
    print(f'Processed {line_count} lines.')

contact_pressure_ = [x for _,x in sorted(zip(x_coord,contact_pressure))]
contact_pressure = contact_pressure_
x_coord.sort()
ax1.plot(x_coord,contact_pressure, 'g--', linewidth=2, label='Dual mortar')

ax1.set_xlabel('Interface longitudinal dimension [$m$]')
ax1.set_ylabel('Contact pressure [$Pa$]')
ax1.grid()
ax1.legend(loc='upper right', fancybox=True, framealpha=0.4)
ax1.set_xlim(-0.075, 0.075)
#ax1.set_ylim(4, 16)

# ax2 = ax1.twinx()
# ax2.set_ylabel('Temperature [$K$]', color='tab:blue')  # we already handled the x-label with ax1


# ax1.show()

# Comparison of stresses
# ax2.plot(x,stress_xy_x_mat, label='Stress from Tensor Mechanics')
#
# ax2.set_xlabel('Z direction')
# ax2.set_ylabel('YZ stress (Pa)')
# ax2.set_title('Material-based periodic solution of bi-phase problem')
# ax2.legend(loc='upper left')
# ax2.grid()

plt.savefig('contact.pdf')
plt.show()
