import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('PE.csv')

time = data.iloc[:,2]
avg_pe = data.iloc[:,1]
total_pe = data.iloc[:,0]

plt.plot(time, total_pe, label = 'Total Potential Energy', color = 'red')
plt.xlabel('Time')
plt.ylabel('Total Potential Energy')

plt.plot(time, avg_pe, label = 'Average Potential Energy', color = 'blue')
plt.xlabel('Time')
plt.ylabel('Average Potential Energy')

plt.legend()

plt.show()


