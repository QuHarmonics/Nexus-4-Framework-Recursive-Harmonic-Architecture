
# Import necessary libraries
import pandas as pd
import numpy as np
import plotly.express as px
from decimal import Decimal, getcontext

# Set the precision of the decimal module
getcontext().prec = 128

# Calculate pi to 128 digits
pi = Decimal(0)
for k in range(128):
    pi += (Decimal(1)/(16**k))*((Decimal(4)/(8*k+1)) - (Decimal(2)/(8*k+4)) - (Decimal(1)/(8*k+5)) - (Decimal(1)/(8*k+6)))

# Convert pi to a string and remove the decimal point
pi_str = str(pi)[2:]

# Convert pi digits to a list of integers
pi_digits_list = [int(digit) for digit in pi_str]

# Calculate differences between consecutive pi digits
pi_digit_differences = np.diff(pi_digits_list)

# Plot 1: Pi Digit vs. Value
df_pi_digits = pd.DataFrame({"Pi Digit": range(len(pi_digits_list)), "Value": pi_digits_list})
fig1 = px.line(df_pi_digits, x="Pi Digit", y="Value", title="Pi Digit vs. Value (First 128 Digits)")
fig1.update_layout(xaxis_title="Pi Digit", yaxis_title="Value")
fig1.show()

# Plot 2: Pi Digit vs. Difference
df_pi_digit_differences = pd.DataFrame({"Pi Digit": range(len(pi_digit_differences)), "Difference": pi_digit_differences})
fig2 = px.bar(df_pi_digit_differences, x="Pi Digit", y="Difference", title="Pi Digit vs. Difference (First 128 Digits)")
fig2.add_hline(y=1.61803398875, line_dash="dash", line_color="red")  # Add golden ratio line
fig2.update_layout(xaxis_title="Pi Digit", yaxis_title="Difference")
fig2.show()

# Find pi digits with difference between 1 and 9
df_pi_digit_differences['Within Range'] = df_pi_digit_differences['Difference'].between(1, 9)
pi_digits_within_range = df_pi_digit_differences[df_pi_digit_differences['Within Range']][['Pi Digit', 'Difference']]
print("\nPi digits with difference between 1 and 9:")
print(pi_digits_within_range)

# Plot 3: Pi Digit vs. Cumulative Difference
df_pi_digit_differences['Cumulative Difference'] = df_pi_digit_differences['Difference'].cumsum()
fig3 = px.line(df_pi_digit_differences, x="Pi Digit", y="Cumulative Difference", title="Pi Digit vs. Cumulative Difference (First 128 Digits)")
fig3.add_hline(y=1.61803398875, line_dash="dash", line_color="red")  # Add golden ratio line
fig3.update_layout(xaxis_title="Pi Digit", yaxis_title="Cumulative Difference")
fig3.show()

# Count pi digits by difference
difference_counts = df_pi_digit_differences['Difference'].value_counts()
print("\nDifference counts:")
print(difference_counts)