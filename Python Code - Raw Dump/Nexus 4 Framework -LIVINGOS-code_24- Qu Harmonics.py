import pandas as pd
import random

# Load symbolic agent learning log generated from Phase V
df_echo = pd.read_csv("echo_agents_symbolic_learning_log.csv")

# Function: simulate symbolic agents predicting ZPHC based on early STI data
def simulate_predictive_agents(df, hops_to_see=5):
    predictions = []

    for i, row in df.iterrows():
        agent_id = row["agent_id"]
        
        # Parse stored trust signature
        trust_signature = eval(row["trust_signature"])  # safely interpret stringified list
        observed = trust_signature[:hops_to_see]
        avg_observed = sum(observed) / len(observed)
        peak = max(observed)
        slope = observed[-1] - observed[0]

        # Simple heuristic prediction model for ZPHC based on early pattern
        zphc_pred = avg_observed >= 0.68 and peak >= 0.72 and slope >= 0

        predictions.append({
            "agent_id": agent_id,
            "packet_id": row["packet_id"],
            "observed_hops": hops_to_see,
            "avg_partial_trust": round(avg_observed, 4),
            "peak_trust": round(peak, 4),
            "slope": round(slope, 4),
            "predicted_zphc": zphc_pred,
            "actual_zphc": row["final_zphc"],
            "correct": zphc_pred == row["final_zphc"]
        })

    return pd.DataFrame(predictions)

# Execute prediction simulation
df_pred = simulate_predictive_agents(df_echo)

# Export result to CSV
df_pred.to_csv("symbolic_predictive_agent_results.csv", index=False)
print("✅ Results saved to: symbolic_predictive_agent_results.csv")
