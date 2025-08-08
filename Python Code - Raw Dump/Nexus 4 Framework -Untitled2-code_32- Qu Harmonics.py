import streamlit as st
import numpy as np

# The long data string
data = (
    "5143614111020411149104143134141356110505813117310111126390279591034141211813121171412512910123456789"
    "0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567894567891011121314150"
    "1234567891011121314150123456789101112131415012345678910111213141501234567891011121314150123456789101112131415"
    # Truncated for brevity, use the full data string
)

# Interactive slider for container width
width = st.slider("Adjust container width:", min_value=10, max_value=100, value=50)

# Split the data string into rows based on the current width
rows = [data[i:i+width] for i in range(0, len(data), width)]

# Display the text in a cascading format
st.write("### Cascading Text:")
st.write("\n".join(rows))
