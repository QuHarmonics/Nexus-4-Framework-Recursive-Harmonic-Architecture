# Prepare summary results for display
summary_df = pd.DataFrame({
    "Iteration": range(1, iterations + 1),
    "Base Lattice Entropy": base_entropy_values,
    "Perturbed Lattice Entropy": perturbed_entropy_values,
})
summary_df["Entropy Difference"] = np.abs(summary_df["Base Lattice Entropy"] - summary_df["Perturbed Lattice Entropy"])

# Display the summary results as a table
print("Recursive Feedback Analysis Summary:")
print(summary_df)

# Optional: Save the summary to a CSV file for further analysis
summary_df.to_csv("recursive_feedback_summary.csv", index=False)

# Plot summary data as a bar chart for visual interpretation
plt.figure(figsize=(12, 6))
plt.bar(summary_df["Iteration"], summary_df["Base Lattice Entropy"], label="Base Lattice Entropy", alpha=0.7)
plt.bar(summary_df["Iteration"], summary_df["Perturbed Lattice Entropy"], label="Perturbed Lattice Entropy", alpha=0.7)
plt.plot(summary_df["Iteration"], summary_df["Entropy Difference"], color="red", marker="o", label="Entropy Difference")
plt.axhline(y=entropy_threshold, color="green", linestyle="--", label="Entropy Threshold")
plt.title("Recursive Feedback Entropy Analysis")
plt.xlabel("Iteration")
plt.ylabel("Entropy")
plt.legend()
plt.grid(True)
plt.show()
