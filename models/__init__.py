"""
CRISPR-Cas9 Simulation Models

This package contains the core models and logic for the CRISPR-Cas9 simulation.
"""

# Import all models and functions to make them available at the package level
# Note: We're using absolute imports to avoid circular imports

# These will be populated by the actual imports below
Cas9 = None
calculate_on_target_score = None
off_target_score = None
pam_score = None
secondary_structure_penalty = None
ChromatinContext = None
cas9_activity = None
repair_pathway_prob = None
simulate_delivery = None
delivery_efficiency = None

# Import everything at the end to avoid circular imports
try:
    from CRISPRcas9_simV3.models.cas9 import Cas9
    from CRISPRcas9_simV3.models.scoring import calculate_on_target_score, off_target_score, pam_score, secondary_structure_penalty
    from CRISPRcas9_simV3.models.simulation import ChromatinContext, cas9_activity, repair_pathway_prob, simulate_delivery, delivery_efficiency
except ImportError:
    # This allows the imports to work when running the package directly
    from .cas9 import Cas9
    from .scoring import calculate_on_target_score, off_target_score, pam_score, secondary_structure_penalty
    from .simulation import ChromatinContext, cas9_activity, repair_pathway_prob, simulate_delivery, delivery_efficiency

__all__ = [
    'Cas9',
    'calculate_on_target_score',
    'off_target_score',
    'pam_score',
    'secondary_structure_penalty',
    'ChromatinContext',
    'cas9_activity',
    'repair_pathway_prob',
    'simulate_delivery',
    'delivery_efficiency'
]
