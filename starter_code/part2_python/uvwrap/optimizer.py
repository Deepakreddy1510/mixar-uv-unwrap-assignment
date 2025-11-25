"""
Parameter optimization using grid search

TEMPLATE - YOU IMPLEMENT

Finds best parameters for a given mesh by testing multiple combinations.
"""

import itertools


def optimize_parameters(mesh_path, target_metric='stretch', verbose=True):
    """
    Find optimal unwrapping parameters for a mesh

    Uses grid search over parameter space to minimize/maximize target metric.

    Args:
        mesh_path: Path to mesh file
        target_metric: Metric to optimize ('stretch', 'coverage', 'angle_distortion')
        verbose: Print progress

    Returns:
        tuple: (best_params, best_score)
            best_params: Dictionary of parameters
            best_score: Best metric value achieved

    IMPLEMENTATION REQUIRED
    """
    # TODO: Implement
    #
    # Algorithm:
    # 1. Define parameter search space:
    #      angle_thresholds = [20, 30, 40, 50]
    #      min_island_sizes = [5, 10, 20, 50]
    #
    # 2. For each parameter combination:
    #      Unwrap mesh with these params
    #      Compute target metric
    #      Track best result
    #
    # 3. Return best parameters and score

    # Parameter search space
    angle_thresholds = [20, 30, 40, 50]  # degrees
    min_island_sizes = [5, 10, 20, 50]   # faces

    best_params = None
    best_score = float('inf') if target_metric != 'coverage' else 0.0

    total_combinations = len(angle_thresholds) * len(min_island_sizes)
    current = 0

    if verbose:
        print(f"Testing {total_combinations} parameter combinations...")
        print(f"Target metric: {target_metric}")
        print()

    # TODO: YOUR CODE HERE
    # Implement grid search

    return best_params, best_score


# Example usage
if __name__ == "__main__":
    # Test optimizer
    # best_params, best_score = optimize_parameters(
    #     "../../test_data/meshes/01_cube.obj",
    #     target_metric='stretch',
    #     verbose=True
    # )
    # print(f"\nBest parameters: {best_params}")
    # print(f"Best score: {best_score:.3f}")
    pass
