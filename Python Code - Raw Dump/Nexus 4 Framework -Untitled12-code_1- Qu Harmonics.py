import numpy as np
import plotly.graph_objects as go

def create_single_matrix(center, frame_size):
    """
    Create a single matrix (box) with a given center and frame size.

    Args:
        center (tuple): The center of the matrix (x, y, z).
        frame_size (int): The size of the matrix.

    Returns:
        dict: Dictionary containing vertices and edges of the matrix.
    """
    half = frame_size / 2
    cx, cy, cz = center

    vertices = np.array([
        [cx - half, cy - half, cz - half],  # Bottom-left-front
        [cx + half, cy - half, cz - half],  # Bottom-right-front
        [cx + half, cy + half, cz - half],  # Top-right-front
        [cx - half, cy + half, cz - half],  # Top-left-front
        [cx - half, cy - half, cz + half],  # Bottom-left-back
        [cx + half, cy - half, cz + half],  # Bottom-right-back
        [cx + half, cy + half, cz + half],  # Top-right-back
        [cx - half, cy + half, cz + half],  # Top-left-back
    ])

    edges = [
        (0, 1), (1, 2), (2, 3), (3, 0),  # Bottom face
        (4, 5), (5, 6), (6, 7), (7, 4),  # Top face
        (0, 4), (1, 5), (2, 6), (3, 7)   # Vertical edges
    ]

    return {'vertices': vertices, 'edges': edges, 'frame_size': frame_size}

def recursive_matrices(frame_size, recursion_depth):
    """
    Generate recursive matrices starting at (0, 0, 0) and using the calculated center of the previous box.

    Args:
        frame_size (int): The base size of the smallest matrix.
        recursion_depth (int): The number of recursive levels.

    Returns:
        list: A list of matrices with their vertices and edges.
    """
    matrices = []
    center = (0, 0, 0)

    for i in range(recursion_depth):
        current_size = frame_size * (2 ** i)
        matrix = create_single_matrix(center, current_size)
        matrices.append(matrix)
        # Calculate the new center as the center of the previous matrix
        center = tuple(np.mean(matrix['vertices'], axis=0))

    return matrices

def visualize_recursive_matrices_plotly(matrices):
    """
    Visualize recursive matrices in 3D using Plotly.

    Args:
        matrices (list): A list of dictionaries containing vertices and edges.
    """
    fig = go.Figure()

    colors = ['red', 'blue', 'green', 'purple', 'orange']

    for i, matrix in enumerate(matrices):
        vertices = matrix['vertices']
        edges = matrix['edges']
        color = colors[i % len(colors)]

        # Add edges as 3D lines
        for edge in edges:
            start, end = edge
            fig.add_trace(go.Scatter3d(
                x=[vertices[start][0], vertices[end][0]],
                y=[vertices[start][1], vertices[end][1]],
                z=[vertices[start][2], vertices[end][2]],
                mode='lines',
                line=dict(color=color, width=3),
                name=f'Box {i + 1}'
            ))

    # Set layout
    fig.update_layout(
        title="Recursive Nested Matrices Aligned Recursively",
        scene=dict(
            xaxis_title="X-axis",
            yaxis_title="Y-axis",
            zaxis_title="Z-axis",
            aspectmode='cube'
        ),
        margin=dict(l=0, r=0, t=50, b=0)
    )

    fig.show()

def main():
    frame_size = 512  # Base size of the smallest matrix
    recursion_depth = 3  # Number of recursive levels

    matrices = recursive_matrices(frame_size, recursion_depth)
    visualize_recursive_matrices_plotly(matrices)

if __name__ == "__main__":
    main()
