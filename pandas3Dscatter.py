import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

df = pd.read_csv('~/a', parse_dates=True)
#print(df.head())
#df['H-L'] = df.High - df.Low
#df['100MA'] = pd.rolling_mean(df['Close'], 100)

threedee = plt.figure().gca(projection='3d')
threedee.scatter(df['x'], df['y'], df['z'])
threedee.set_xlabel('x')
threedee.set_ylabel('y')
threedee.set_zlabel('z')
plt.show()
