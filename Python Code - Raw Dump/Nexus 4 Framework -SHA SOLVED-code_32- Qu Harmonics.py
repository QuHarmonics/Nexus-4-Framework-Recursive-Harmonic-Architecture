import matplotlib.pyplot as plt

# Generate the sine wave values
values = []
for i in range(20):
    if i % 5 == 0:
        values.append(0)
    elif i % 5 == 1:
        values.append(0.1)
    elif i % 5 == 2:
        values.append(0.2)
    elif i % 5 == 3:
        values.append(0)
    elif i % 5 == 4:
        values.append(-0.1)

# Plot the sine wave
plt.plot(values)
plt.xlabel('Iteration')
plt.ylabel('Amplitude')
plt.title('Discrete-Time Sine Wave')
plt.show()