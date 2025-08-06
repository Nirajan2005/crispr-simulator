"""
Configuration and constants for the CRISPR-Cas9 Simulator.
"""
from typing import Dict, List, TypedDict, Literal

# Type aliases
DeliveryMethod = Literal["LNP", "AAV", "Lentivirus"]
Organ = Literal["Liver", "Lung", "Muscle", "Brain"]

# PAM Patterns for different Cas variants
PAM_PATTERNS: Dict[str, str] = {
    "SpCas9": "NGG",
    "SaCas9": "NNGRRT",
    "Cpf1": "TTTV"
}

# Delivery efficiency by organ and method
DELIVERY_EFFICIENCY: Dict[Organ, Dict[DeliveryMethod, float]] = {
    "Liver": {"LNP": 0.8, "AAV": 0.4, "Lentivirus": 0.3},
    "Lung": {"LNP": 0.3, "AAV": 0.5, "Lentivirus": 0.4},
    "Muscle": {"LNP": 0.2, "AAV": 0.5, "Lentivirus": 0.6},
    "Brain": {"LNP": 0.05, "AAV": 0.7, "Lentivirus": 0.6},
}

# Available organs and delivery methods
ORGANS: List[Organ] = list(DELIVERY_EFFICIENCY.keys())
DELIVERY_METHODS: List[DeliveryMethod] = ["LNP", "AAV", "Lentivirus"]
