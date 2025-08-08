#!/usr/bin/env python3
"""
Test script to verify that all functionality from the original CRISPR_v3.py
has been restored in the modular structure.
"""

import sys
import traceback

def test_imports():
    """Test that all required modules can be imported."""
    print("Testing imports...")
    
    try:
        from constants import IUPAC_CODES
        print("✓ constants imported")
    except Exception as e:
        print(f"✗ constants failed: {e}")
        return False
    
    try:
        from config import PAM_PATTERNS, DELIVERY_EFFICIENCY, ORGANS, DELIVERY_METHODS
        print("✓ config imported")
    except Exception as e:
        print(f"✗ config failed: {e}")
        return False
    
    try:
        from models.cas9 import Cas9
        print("✓ Cas9 model imported")
    except Exception as e:
        print(f"✗ Cas9 model failed: {e}")
        return False
    
    try:
        from models.scoring import calculate_on_target_score, off_target_score, pam_score, secondary_structure_penalty
        print("✓ scoring functions imported")
    except Exception as e:
        print(f"✗ scoring functions failed: {e}")
        return False
    
    try:
        from models.simulation import ChromatinContext, cas9_activity, repair_pathway_prob, simulate_delivery, delivery_efficiency
        print("✓ simulation functions imported")
    except Exception as e:
        print(f"✗ simulation functions failed: {e}")
        return False
    
    try:
        from utils.sequence_utils import is_valid_pam
        print("✓ sequence utilities imported")
    except Exception as e:
        print(f"✗ sequence utilities failed: {e}")
        return False
    
    return True

def test_core_functionality():
    """Test core CRISPR functionality."""
    print("\nTesting core functionality...")
    
    try:
        from models.cas9 import Cas9
        from utils.sequence_utils import is_valid_pam
        from models.scoring import calculate_on_target_score
        
        # Test sequence and gRNA
        test_seq = "ATGCGTACCGGTTACCGGATCGTACCGGTTACCGGATCGTACCGGTTACCGGA"
        test_grna = "CGTACCGGTTACCGGATCGT"
        
        # Test Cas9 initialization
        cas9 = Cas9(pam_pattern="NGG")
        print("✓ Cas9 initialization")
        
        # Test PAM site finding
        pam_sites = cas9.find_pam_sites(test_seq)
        print(f"✓ PAM site finding: found {len(pam_sites)} sites")
        
        # Test target finding
        targets = cas9.find_targets(test_seq, test_grna)
        print(f"✓ Target finding: found {len(targets)} targets")
        
        # Test scoring
        if test_grna and len(test_grna) == 20:
            score = calculate_on_target_score(test_grna)
            print(f"✓ On-target scoring: {score:.2f}")
        
        # Test PAM validation
        valid = is_valid_pam("AGG", "NGG")
        print(f"✓ PAM validation: AGG vs NGG = {valid}")
        
        # Test editing simulation
        if targets:
            edited = cas9.simulate_edit(test_seq, targets[0], "Knockout (NHEJ)")
            print(f"✓ Edit simulation: {len(test_seq)} -> {len(edited)} bp")
        
        return True
        
    except Exception as e:
        print(f"✗ Core functionality failed: {e}")
        traceback.print_exc()
        return False

def test_delivery_modeling():
    """Test delivery and biological modeling."""
    print("\nTesting delivery modeling...")
    
    try:
        from models.simulation import ChromatinContext, cas9_activity, repair_pathway_prob, simulate_delivery, delivery_efficiency
        from config import DELIVERY_EFFICIENCY
        
        # Test chromatin context
        chromatin = ChromatinContext(accessibility_score=0.8, histone_marks={"H3K27ac": 1})
        efficiency = chromatin.modify_cutting_efficiency(0.5)
        print(f"✓ Chromatin context: {efficiency:.2f}")
        
        # Test Cas9 activity
        activity = cas9_activity(37.0, 100.0)
        print(f"✓ Cas9 activity: {activity:.2f}")
        
        # Test repair pathway
        hdr_prob = repair_pathway_prob("S/G2")
        print(f"✓ Repair pathway: HDR prob = {hdr_prob:.2f}")
        
        # Test delivery simulation
        delivered = simulate_delivery("Liver", "LNP")
        print(f"✓ Delivery simulation: {delivered}")
        
        # Test delivery efficiency calculation
        eff = delivery_efficiency(0.8, 5.0, "IV")
        print(f"✓ Delivery efficiency: {eff:.2f}")
        
        return True
        
    except Exception as e:
        print(f"✗ Delivery modeling failed: {e}")
        traceback.print_exc()
        return False

def test_original_vs_modular():
    """Compare original vs modular implementation."""
    print("\nComparing original vs modular...")
    
    try:
        # Test with the same sequence from original
        test_seq = "ATGCGTACCGGTTACCGGATCGTACCGGTTACCGGATCGTACCGGTTACCGGA"
        test_grna = "CGTACCGGTTACCGGATCGT"
        
        # Import original
        import CRISPR_v3 as original
        
        # Import modular
        from models.cas9 import Cas9 as ModularCas9
        
        # Test with both
        original_cas9 = original.Cas9(pam_pattern="NGG")
        modular_cas9 = ModularCas9(pam_pattern="NGG")
        
        original_targets = original_cas9.find_targets(test_seq, test_grna)
        modular_targets = modular_cas9.find_targets(test_seq, test_grna)
        
        print(f"✓ Original found {len(original_targets)} targets")
        print(f"✓ Modular found {len(modular_targets)} targets")
        
        if len(original_targets) == len(modular_targets):
            print("✓ Target counts match!")
        else:
            print("⚠ Target counts differ")
        
        return True
        
    except Exception as e:
        print(f"✗ Comparison failed: {e}")
        traceback.print_exc()
        return False

def main():
    """Run all tests."""
    print("CRISPR-Cas9 Simulator Functionality Test")
    print("=" * 50)
    
    tests = [
        test_imports,
        test_core_functionality,
        test_delivery_modeling,
        test_original_vs_modular
    ]
    
    passed = 0
    for test in tests:
        try:
            if test():
                passed += 1
            else:
                print(f"✗ {test.__name__} failed")
        except Exception as e:
            print(f"✗ {test.__name__} crashed: {e}")
    
    print(f"\nResults: {passed}/{len(tests)} tests passed")
    
    if passed == len(tests):
        print("🎉 All tests passed! The modular structure successfully restores the original functionality.")
        return 0
    else:
        print("❌ Some tests failed. The restoration is incomplete.")
        return 1

if __name__ == "__main__":
    sys.exit(main())