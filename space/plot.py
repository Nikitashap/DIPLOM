import numpy as np
import matplotlib.pyplot as plt

# Чтение данных из файла ress.txt
data = np.loadtxt('res2.txt')

# Разделяем данные на t, y_explicit, y_implicit, y_exact
t = data[:, 0]
x = data[:, 1]
y_implicit = data[:, 2]
y_exact = data[:, 3]


# Построение графиков
plt.figure(figsize=(10, 6))

# # Построение для явной схемы
# plt.plot(t, y_explicit, label='Explicit Scheme', color='blue')

# # Построение для неявной схемы

x1 = np.linspace(0,1,100)
y1 = np.linspace(2.71828,0,100) 
for i in range(0,100,10):
	print(abs(y_exact[i]-y1[i]), abs(y_implicit[i]-y1[i]))

plt.plot(x, y_implicit, label='Implicit Scheme', color='green')

y2=[abs(y_exact[i]-y_implicit[i]) for i in range(len(y_exact))]

# Построение для точного решения
plt.plot(x, y_exact, label='Exact Solution', linestyle='--', color='red')
# plt.plot(x1,y1)



# Настройка графика
plt.xlabel('x')
plt.ylabel('u(x,1)')
plt.title('Numerical and Exact Solutions')
plt.legend()
plt.grid(True)

# Показываем график
plt.show()

plt.plot(x,y2)

plt.xlabel('Time (t)')
plt.ylabel('y(t)')
plt.title('Error')
plt.legend()
plt.grid(True)

# Показываем график
plt.show()

