import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

# Extracting the instructions from the provided data
data = """
41 inc ecx
54 push esp
47 inc edi
43 inc ebx
"""

# Parsing the data into a structured format
lines = data.strip().split("\n")
instructions = [line.split()[1] for line in lines]

# Counting occurrences of each instruction
instruction_counts = Counter(instructions)

# Creating a DataFrame for visualization and further analysis
instruction_df = pd.DataFrame.from_dict(instruction_counts, orient='index', columns=['Count'])
instruction_df.index.name = 'Instruction'
instruction_df = instruction_df.sort_values(by='Count', ascending=False)

# Displaying the data to the user
import ace_tools as tools; tools.display_dataframe_to_user(name="Instruction Count Map", dataframe=instruction_df)

# Plotting the instruction frequency
plt.figure(figsize=(10, 6))
instruction_df.plot(kind='bar', legend=False, figsize=(10, 6))
plt.title("Instruction Frequency Distribution")
plt.ylabel("Count")
plt.xlabel("Instruction")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()
