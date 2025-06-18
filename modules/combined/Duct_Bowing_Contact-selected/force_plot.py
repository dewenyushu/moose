import pandas as pd
import matplotlib.pyplot as plt

default_path = "./default_8proc_outputs/tlp_force_plots_finer_lps.csv"
default_df  = pd.read_csv(default_path, sep=',')

pp_change_path = "./ppchange_8proc_outputs/tlp_force_plots_finer_lps.csv"
pp_change_df = pd.read_csv(pp_change_path, sep=',')

column_names = pp_change_df.columns


for col in column_names:
  plt.figure()
  if (col!='time'):
    plt.plot(pp_change_df['time'], pp_change_df[col], '-', label = col+"_pp")
    plt.plot(default_df['time'], default_df[col], '--', label = col)

  # Adding labels and title
  plt.xlabel('Time')
  plt.ylabel('Contact Force')
  plt.xlim([0, 1])
  # plt.ylim([0,800])
  # plt.title('Contact Force Compare')
  plt.legend()
  plt.grid(True)
  # plt.xscale('log', base=2)
  # plt.yscale('log', base=10)

  # Save the plot to a PNG file
  plt.savefig('Figures/contact_force_'+col+'.png')
