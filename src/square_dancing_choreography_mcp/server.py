"""
Square Dancing MCP Server - Extended Categorical Edition

Adds three major capabilities:
1. Generative Effects Analysis - Measure compositional emergence
2. Hip Hop Cascade Integration - Receive rhythmic parameters from hip hop server
3. Comparative Visualization - Demonstrate isomorphism with aesthetic composition

Core insight: Square dancing is physically grounded proof that categorical
composition produces genuine generative effects.
"""

from fastmcp import FastMCP
from typing import Dict, List, Optional, Literal, Tuple
import yaml
import math
from pathlib import Path
import json

# Initialize FastMCP server
mcp = FastMCP("squaredance-choreography-extended")

# Load YAML olog specification (using embedded version for now)
OLOG_YAML = """
olog_name: squaredance-categorical
version: "2.0"

dimensions:
  step_pattern:
    semantic_meaning: "Basic movement vocabulary"
    tokens: [walking_step, shuffle_step, quick_allemande, slow_promenade,
             balance_step, swing_step, do_si_do, grand_right_left, complex_figure]
    
  formation_type:
    semantic_meaning: "Spatial arrangement of dancers"
    tokens: [square_set, circle_formation, longways_set, scatter_promenade,
             triangle_set, couple_facing, star_formation, lines_of_four, 
             ocean_wave, facing_lines, right_left_grand_circle]
    
  timing_structure:
    semantic_meaning: "Temporal organization of calls"
    tokens: [sixteen_beat_figure, thirty_two_beat_chorus, sixty_four_beat_sequence,
             one_twenty_eight_beat_extended, variable_patter, singing_call_structure]
    
  call_style:
    semantic_meaning: "Energy and tempo of calling"
    tokens: [smooth_promenade, moderate_pace, accelerating_figures,
             swing_your_partner, promenade_home, syncopated_hash, flowing_continuous]
    
  couple_dynamics:
    semantic_meaning: "How couples move through formations"
    tokens: [lead_couple_position, visiting_couples, all_eight_circulate,
             heads_sides_alternating, corner_progression, partner_swing,
             courtesy_turn, allemande_rotation]
    
  caller_interaction:
    semantic_meaning: "Style of caller direction"
    tokens: [standard_calls, patter_calling, singing_calls, 
             competitive_hash, cooperative_flow, workshop_teaching]
    
  floor_participation:
    semantic_meaning: "How many squares/dancers involved"
    tokens: [single_square, demonstration_squares, all_squares_active,
             progressive_squares, scattered_minglers]
    
  geometric_symmetry:
    semantic_meaning: "Symmetry of formation geometry"
    tokens: [four_fold_square, eight_fold_circle, bilateral_lines,
             radial_star, translational_chain, point_symmetric, asymmetric]
    
  call_complexity:
    semantic_meaning: "Difficulty level of figures"
    tokens: [basic_mainstream, plus_level, advanced_challenging,
             c1_complex, experimental_choreography]
             
  flow_quality:
    semantic_meaning: "Emergent smoothness of call transitions"
    tokens: [choppy_disjointed, adequate_functional, smooth_flowing, 
             seamless_elegant, virtuosic_effortless]

# Call Library - Minimal but complete enough for demonstrations
call_library:
  # ---- BASIC / MAINSTREAM ----
  allemande_left:
    timing: 8
    movement_type: turn
    preconditions: [square_set, ocean_wave]
    postconditions: [ocean_wave, right_left_grand_circle]
    symmetry_effect: maintains
    geometric_primitive: arc
    level: basic
    
  do_si_do:
    timing: 8
    movement_type: walk
    preconditions: [square_set, couple_facing]
    postconditions: [square_set, couple_facing]
    symmetry_effect: maintains
    geometric_primitive: circle
    level: basic
    
  swing_thru:
    timing: 8
    movement_type: weave
    preconditions: [ocean_wave]
    postconditions: [ocean_wave]
    symmetry_effect: transforms
    geometric_primitive: weave
    level: basic

  promenade:
    timing: 16
    movement_type: walk
    preconditions: [square_set, couple_facing]
    postconditions: [square_set]
    symmetry_effect: maintains
    geometric_primitive: circle
    level: basic
    
  square_thru:
    timing: 16
    movement_type: weave
    preconditions: [facing_lines, square_set]
    postconditions: [facing_lines]
    symmetry_effect: breaks
    geometric_primitive: weave
    level: basic
    
  circle_left:
    timing: 8
    movement_type: circle
    preconditions: [square_set]
    postconditions: [circle_formation]
    symmetry_effect: transforms
    geometric_primitive: circle
    level: basic
    
  circle_right:
    timing: 8
    movement_type: circle
    preconditions: [square_set]
    postconditions: [circle_formation]
    symmetry_effect: transforms
    geometric_primitive: circle
    level: basic
    
  pass_thru:
    timing: 4
    movement_type: walk
    preconditions: [facing_lines, couple_facing]
    postconditions: [facing_lines]
    symmetry_effect: breaks
    geometric_primitive: line
    level: basic
    
  star_right:
    timing: 8
    movement_type: rotation
    preconditions: [square_set, circle_formation]
    postconditions: [star_formation]
    symmetry_effect: transforms
    geometric_primitive: star
    level: basic

  star_left:
    timing: 8
    movement_type: rotation
    preconditions: [square_set, circle_formation]
    postconditions: [star_formation]
    symmetry_effect: transforms
    geometric_primitive: star
    level: basic

  right_and_left_thru:
    timing: 8
    movement_type: turn
    preconditions: [facing_lines, couple_facing]
    postconditions: [facing_lines, couple_facing]
    symmetry_effect: maintains
    geometric_primitive: line
    level: basic

  ladies_chain:
    timing: 8
    movement_type: turn
    preconditions: [facing_lines, square_set]
    postconditions: [facing_lines, square_set]
    symmetry_effect: maintains
    geometric_primitive: arc
    level: basic

  courtesy_turn:
    timing: 4
    movement_type: turn
    preconditions: [couple_facing]
    postconditions: [couple_facing]
    symmetry_effect: maintains
    geometric_primitive: arc
    level: basic

  california_twirl:
    timing: 4
    movement_type: turn
    preconditions: [couple_facing]
    postconditions: [couple_facing]
    symmetry_effect: maintains
    geometric_primitive: arc
    level: basic

  dive_thru:
    timing: 4
    movement_type: walk
    preconditions: [facing_lines]
    postconditions: [square_set]
    symmetry_effect: restores
    geometric_primitive: line
    level: basic

  bend_the_line:
    timing: 4
    movement_type: walk
    preconditions: [facing_lines, lines_of_four]
    postconditions: [couple_facing]
    symmetry_effect: transforms
    geometric_primitive: line
    level: basic

  slide_thru:
    timing: 4
    movement_type: walk
    preconditions: [couple_facing, facing_lines]
    postconditions: [couple_facing]
    symmetry_effect: maintains
    geometric_primitive: line
    level: basic

  grand_right_and_left:
    timing: 16
    movement_type: weave
    preconditions: [right_left_grand_circle]
    postconditions: [square_set]
    symmetry_effect: restores
    geometric_primitive: circle
    level: basic

  star_promenade:
    timing: 8
    movement_type: walk
    preconditions: [star_formation]
    postconditions: [square_set]
    symmetry_effect: restores
    geometric_primitive: star
    level: basic

  backtrack:
    timing: 4
    movement_type: walk
    preconditions: [star_formation]
    postconditions: [circle_formation]
    symmetry_effect: transforms
    geometric_primitive: circle
    level: basic

  # ---- PLUS LEVEL ----
  touch_a_quarter:
    timing: 4
    movement_type: turn
    preconditions: [couple_facing, facing_lines]
    postconditions: [ocean_wave]
    symmetry_effect: transforms
    geometric_primitive: arc
    level: plus

  spin_the_top:
    timing: 8
    movement_type: weave
    preconditions: [ocean_wave]
    postconditions: [ocean_wave]
    symmetry_effect: transforms
    geometric_primitive: weave
    level: plus

  trade:
    timing: 4
    movement_type: turn
    preconditions: [ocean_wave, facing_lines]
    postconditions: [ocean_wave, facing_lines]
    symmetry_effect: maintains
    geometric_primitive: arc
    level: plus

  run:
    timing: 4
    movement_type: walk
    preconditions: [ocean_wave]
    postconditions: [couple_facing]
    symmetry_effect: breaks
    geometric_primitive: arc
    level: plus

  circulate:
    timing: 4
    movement_type: walk
    preconditions: [ocean_wave, lines_of_four]
    postconditions: [ocean_wave, lines_of_four]
    symmetry_effect: maintains
    geometric_primitive: line
    level: plus

  recycle:
    timing: 4
    movement_type: walk
    preconditions: [ocean_wave]
    postconditions: [couple_facing]
    symmetry_effect: breaks
    geometric_primitive: arc
    level: plus

  ferris_wheel:
    timing: 8
    movement_type: circle
    preconditions: [facing_lines]
    postconditions: [couple_facing]
    symmetry_effect: transforms
    geometric_primitive: circle
    level: plus

  flutter_wheel:
    timing: 8
    movement_type: rotation
    preconditions: [couple_facing, facing_lines]
    postconditions: [couple_facing]
    symmetry_effect: maintains
    geometric_primitive: circle
    level: plus

  sweep_a_quarter:
    timing: 4
    movement_type: walk
    preconditions: [couple_facing]
    postconditions: [facing_lines]
    symmetry_effect: transforms
    geometric_primitive: arc
    level: plus

  scoot_back:
    timing: 8
    movement_type: walk
    preconditions: [ocean_wave]
    postconditions: [ocean_wave]
    symmetry_effect: maintains
    geometric_primitive: line
    level: plus

  walk_and_dodge:
    timing: 4
    movement_type: walk
    preconditions: [ocean_wave]
    postconditions: [couple_facing]
    symmetry_effect: breaks
    geometric_primitive: line
    level: plus

  zoom:
    timing: 4
    movement_type: walk
    preconditions: [lines_of_four, facing_lines]
    postconditions: [lines_of_four]
    symmetry_effect: maintains
    geometric_primitive: line
    level: plus

  # ---- ADVANCED / CHALLENGE ----
  spin_chain_thru:
    timing: 16
    movement_type: weave
    preconditions: [ocean_wave]
    postconditions: [ocean_wave]
    symmetry_effect: transforms
    geometric_primitive: weave
    level: advanced

  relay_the_deucey:
    timing: 16
    movement_type: weave
    preconditions: [ocean_wave]
    postconditions: [ocean_wave]
    symmetry_effect: transforms
    geometric_primitive: weave
    level: advanced

  load_the_boat:
    timing: 16
    movement_type: weave
    preconditions: [facing_lines]
    postconditions: [facing_lines]
    symmetry_effect: breaks
    geometric_primitive: weave
    level: advanced

  eight_chain_thru:
    timing: 16
    movement_type: weave
    preconditions: [facing_lines]
    postconditions: [facing_lines]
    symmetry_effect: breaks
    geometric_primitive: weave
    level: advanced

  extend:
    timing: 4
    movement_type: walk
    preconditions: [ocean_wave]
    postconditions: [facing_lines]
    symmetry_effect: transforms
    geometric_primitive: line
    level: advanced

  hinge:
    timing: 4
    movement_type: turn
    preconditions: [ocean_wave, lines_of_four]
    postconditions: [ocean_wave]
    symmetry_effect: maintains
    geometric_primitive: arc
    level: advanced

  diamond_circulate:
    timing: 4
    movement_type: walk
    preconditions: [star_formation]
    postconditions: [star_formation]
    symmetry_effect: maintains
    geometric_primitive: star
    level: advanced

  flip_the_diamond:
    timing: 4
    movement_type: walk
    preconditions: [star_formation]
    postconditions: [ocean_wave]
    symmetry_effect: transforms
    geometric_primitive: star
    level: advanced

  cut_the_diamond:
    timing: 4
    movement_type: walk
    preconditions: [star_formation]
    postconditions: [facing_lines]
    symmetry_effect: breaks
    geometric_primitive: star
    level: advanced

  peel_off:
    timing: 4
    movement_type: walk
    preconditions: [lines_of_four, facing_lines]
    postconditions: [lines_of_four]
    symmetry_effect: transforms
    geometric_primitive: line
    level: advanced

  transfer_the_column:
    timing: 16
    movement_type: weave
    preconditions: [lines_of_four]
    postconditions: [ocean_wave]
    symmetry_effect: transforms
    geometric_primitive: weave
    level: advanced

categorical_structures:
  promenade_monoid:
    description: "Basic promenade: identity (stand) + walk cycle + associativity"
    category_type: monoid
    form: "Walk 4 steps, repeat - any grouping of steps works"
    maps_to:
      step_pattern: walking_step
      formation_type: square_set
      timing_structure: sixteen_beat_figure
      call_style: smooth_promenade
      couple_dynamics: partner_swing
      caller_interaction: standard_calls
      floor_participation: single_square
      geometric_symmetry: four_fold_square
      call_complexity: basic_mainstream
    
  square_set_groupoid:
    description: "Square formation with interchangeable couple positions"
    category_type: groupoid
    form: "4 couples (8 dancers), all positions equivalent, invertible rotations"
    maps_to:
      step_pattern: shuffle_step
      formation_type: square_set
      timing_structure: sixty_four_beat_sequence
      call_style: moderate_pace
      couple_dynamics: visiting_couples
      caller_interaction: cooperative_flow
      floor_participation: demonstration_squares
      geometric_symmetry: four_fold_square
      call_complexity: plus_level
    
  sixty_four_beat_sequence:
    description: "Complete square dance sequence with phrase structure"
    category_type: category
    form: "Opener (16) → Figure (32) → Break (8) → Closer (8) = 64 beats"
    maps_to:
      step_pattern: complex_figure
      formation_type: square_set
      timing_structure: sixty_four_beat_sequence
      call_style: flowing_continuous
      couple_dynamics: all_eight_circulate
      caller_interaction: patter_calling
      floor_participation: all_squares_active
      geometric_symmetry: four_fold_square
      call_complexity: advanced_challenging
    
  hash_to_flow_functor:
    description: "Transform competitive hash to cooperative flow preserving structure"
    category_type: functor
    form: "F: CompetitiveHash → CooperativeFlow preserving calls, timing, formations"
    maps_to:
      step_pattern: do_si_do
      formation_type: circle_formation
      timing_structure: variable_patter
      call_style: syncopated_hash
      couple_dynamics: corner_progression
      caller_interaction: competitive_hash
      floor_participation: scattered_minglers
      geometric_symmetry: eight_fold_circle
      call_complexity: c1_complex
"""

# Parse YAML
OLOG_SPEC = yaml.safe_load(OLOG_YAML)
DIMENSIONS = OLOG_SPEC['dimensions']
CALL_LIBRARY = OLOG_SPEC['call_library']
CATEGORICAL_STRUCTURES = OLOG_SPEC['categorical_structures']

# ============================================================================
# ORIGINAL TOOLS (from base server)
# ============================================================================

@mcp.tool()
def get_categorical_structure(structure_name: str) -> Dict:
    """Get categorical structure specification."""
    if structure_name not in CATEGORICAL_STRUCTURES:
        raise ValueError(f"Unknown structure: {structure_name}")
    return CATEGORICAL_STRUCTURES[structure_name]

@mcp.tool()
def build_choreography_description(params: Dict[str, str], structure_name: Optional[str] = None) -> str:
    """Build detailed choreography description from parameters."""
    parts = []
    
    formation = params.get('formation_type', 'square_set')
    parts.append(f"{formation.replace('_', ' ')}")
    
    symmetry = params.get('geometric_symmetry', 'four_fold_square')
    parts.append(f"with {symmetry.replace('_', ' ')}")
    
    step = params.get('step_pattern', 'walking_step')
    parts.append(f"{step.replace('_', ' ')} movements")
    
    timing = params.get('timing_structure', 'sixty_four_beat_sequence')
    parts.append(f"{timing.replace('_', ' ')} timing")
    
    dynamics = params.get('couple_dynamics', 'visiting_couples')
    parts.append(f"{dynamics.replace('_', ' ')} pattern")
    
    count = params.get('formation_count', 8)
    parts.append(f"{count} dancers total")
    
    call = params.get('call_style', 'moderate_pace')
    parts.append(f"{call.replace('_', ' ')} caller style")
    
    description = ", ".join(parts)
    
    if structure_name and structure_name in CATEGORICAL_STRUCTURES:
        structure = CATEGORICAL_STRUCTURES[structure_name]
        description = f"{structure['form']} - {description}"
    
    return description

@mcp.tool()
def build_visualization_prompt(
    params: Dict[str, str],
    structure_name: Optional[str] = None,
    include_notation: bool = True
) -> str:
    """Build complete visualization prompt for choreography diagram."""
    choreo_desc = build_choreography_description(params, structure_name)
    
    viz_elements = [
        "square dance choreography diagram",
        "top-down formation view"
    ]
    
    formation = params.get('formation_type', 'square_set')
    if 'square' in formation:
        viz_elements.append("4 couples in square geometry")
        viz_elements.append("corner positions marked")
    elif 'circle' in formation:
        viz_elements.append("dancers in circle formation")
        viz_elements.append("radial symmetry visible")
    elif 'line' in formation:
        viz_elements.append("dancers in parallel lines")
        viz_elements.append("bilateral symmetry")
    
    symmetry = params.get('geometric_symmetry', 'four_fold_square')
    if 'four_fold' in symmetry:
        viz_elements.append("4-fold rotational symmetry lines")
    elif 'eight_fold' in symmetry:
        viz_elements.append("8-fold radial symmetry")
    elif 'bilateral' in symmetry:
        viz_elements.append("bilateral symmetry axis")
    
    dynamics = params.get('couple_dynamics', 'visiting_couples')
    viz_elements.append(f"movement arrows showing {dynamics.replace('_', ' ')}")
    
    if include_notation:
        timing = params.get('timing_structure', 'sixty_four_beat_sequence')
        viz_elements.append(f"beat markers for {timing.replace('_', ' ')}")
        viz_elements.append("caller notation labels")
        viz_elements.append("phrase boundary indicators")
    
    viz_elements.append("couple numbers (1-4)")
    viz_elements.append("position labels (heads/sides)")
    
    call = params.get('call_style', 'moderate_pace')
    viz_elements.append(f"{call.replace('_', ' ')} tempo indication")
    
    prompt = ", ".join(viz_elements)
    prompt_with_context = f"{prompt}, representing {choreo_desc}"
    
    return prompt_with_context

@mcp.tool()
def get_dimensions() -> Dict:
    """Get all available dimensions and their tokens."""
    return {dim: data['tokens'] for dim, data in DIMENSIONS.items()}

@mcp.tool()
def get_structures() -> Dict:
    """Get available categorical structures."""
    return {
        name: {
            "description": info['description'],
            "category_type": info['category_type'],
            "form": info['form']
        }
        for name, info in CATEGORICAL_STRUCTURES.items()
    }

# ============================================================================
# PHASE 1: GENERATIVE EFFECTS ANALYSIS
# ============================================================================

@mcp.tool()
def detect_generative_effects(
    call_sequence: List[str],
    starting_formation: str = "square_set"
) -> Dict:
    """
    Execute a call sequence and identify emergent properties.
    
    This implements the Φ(s1⊕s2) ≠ Φ(s1)∨Φ(s2) formula from the paper.
    
    Returns:
        - individual_properties: Properties each call would have in isolation
        - composed_properties: Properties of the full sequence
        - emergent_properties: Properties that emerged from composition
        - generative_score: 0.0-1.0 measure of emergence strength
    """
    # Validate calls
    for call in call_sequence:
        if call not in CALL_LIBRARY:
            return {"error": f"Unknown call: {call}"}
    
    # Individual call properties
    individual_properties = []
    for call in call_sequence:
        call_data = CALL_LIBRARY[call]
        individual_properties.append({
            "call": call,
            "timing": call_data['timing'],
            "movement_type": call_data['movement_type'],
            "symmetry_effect": call_data['symmetry_effect'],
            "geometric_primitive": call_data['geometric_primitive']
        })
    
    # Simulate composition through formation trace
    current_formation = starting_formation
    formation_trace = [starting_formation]
    symmetry_trace = [_infer_symmetry(starting_formation)]
    valid_sequence = True
    breaks = []
    
    for i, call in enumerate(call_sequence):
        call_data = CALL_LIBRARY[call]
        
        # Check if call is valid from current formation
        if current_formation not in call_data['preconditions']:
            valid_sequence = False
            breaks.append({
                "at_call": i,
                "call_name": call,
                "issue": f"Cannot execute {call} from {current_formation}",
                "required": call_data['preconditions']
            })
            # Continue simulation anyway for analysis
        
        # Update formation (take first postcondition)
        current_formation = call_data['postconditions'][0]
        formation_trace.append(current_formation)
        
        # Track symmetry evolution
        new_symmetry = _apply_symmetry_effect(
            symmetry_trace[-1],
            call_data['symmetry_effect']
        )
        symmetry_trace.append(new_symmetry)
    
    # Calculate composed properties (emergent from whole sequence)
    composed_properties = {
        "total_timing": sum(CALL_LIBRARY[c]['timing'] for c in call_sequence),
        "formation_trace": formation_trace,
        "symmetry_narrative": symmetry_trace,
        "symmetry_journey": _describe_symmetry_journey(symmetry_trace),
        "flow_quality": _calculate_flow_quality(call_sequence, formation_trace),
        "spatial_coherence": _calculate_spatial_coherence(formation_trace),
        "valid_sequence": valid_sequence,
        "constraint_violations": breaks
    }
    
    # Identify emergent properties
    emergent_properties = []
    
    # 1. Flow quality is purely emergent (not property of individual calls)
    emergent_properties.append({
        "property": "flow_quality",
        "value": composed_properties["flow_quality"],
        "explanation": "Smoothness of transitions emerges from call composition, not from individual calls"
    })
    
    # 2. Symmetry narrative arc
    if len(set(symmetry_trace)) > 1:  # Symmetry changes
        emergent_properties.append({
            "property": "symmetry_arc",
            "value": composed_properties["symmetry_journey"],
            "explanation": "Symmetry breaking and restoration is compositional property"
        })
    
    # 3. Spatial return or departure
    if formation_trace[0] == formation_trace[-1]:
        emergent_properties.append({
            "property": "spatial_return",
            "value": True,
            "explanation": "Sequence returns to starting formation - emergent cycle property"
        })
    else:
        emergent_properties.append({
            "property": "spatial_transformation",
            "value": f"{formation_trace[0]} → {formation_trace[-1]}",
            "explanation": "Net spatial displacement emerges from composition"
        })
    
    # 4. Complexity interaction
    movement_types = [CALL_LIBRARY[c]['movement_type'] for c in call_sequence]
    if len(set(movement_types)) > 1:
        emergent_properties.append({
            "property": "movement_variety",
            "value": list(set(movement_types)),
            "explanation": "Choreographic richness from combining different movement vocabularies"
        })
    
    # Calculate generative score
    # Score based on: flow quality + symmetry changes + emergent properties
    score_components = []
    
    # Flow quality contribution (0-0.3)
    flow_scores = {"choppy": 0.0, "adequate": 0.1, "smooth": 0.2, "seamless": 0.3}
    flow_val = composed_properties["flow_quality"].split("_")[0]
    score_components.append(flow_scores.get(flow_val, 0.15))
    
    # Symmetry narrative contribution (0-0.3)
    symmetry_changes = len(set(symmetry_trace)) - 1
    score_components.append(min(symmetry_changes * 0.1, 0.3))
    
    # Emergent properties count contribution (0-0.4)
    score_components.append(min(len(emergent_properties) * 0.1, 0.4))
    
    generative_score = sum(score_components)
    
    return {
        "call_sequence": call_sequence,
        "starting_formation": starting_formation,
        "individual_properties": individual_properties,
        "composed_properties": composed_properties,
        "emergent_properties": emergent_properties,
        "generative_score": round(generative_score, 2),
        "interpretation": _interpret_generative_score(generative_score)
    }

@mcp.tool()
def analyze_flow_quality(
    call_sequence: List[str],
    starting_formation: str = "square_set"
) -> Dict:
    """
    Measure the "flow quality" that emerges from call composition.
    
    This is the square dance equivalent of convergent_features in tomographic integration!
    
    Returns:
        - smoothness_score: How fluidly calls connect (0-1)
        - symmetry_preservation: Does sequence maintain/break/restore symmetry?
        - spatial_coherence: Formation integrity through sequence (0-1)
        - temporal_coherence: Beat alignment and phrase structure
        - flow_classification: One of [choppy, adequate, smooth, seamless, virtuosic]
    """
    # Validate sequence
    for call in call_sequence:
        if call not in CALL_LIBRARY:
            return {"error": f"Unknown call: {call}"}
    
    # Simulate execution
    current_formation = starting_formation
    formation_trace = [starting_formation]
    transition_scores = []
    
    for call in call_sequence:
        call_data = CALL_LIBRARY[call]
        
        # Score this transition
        if current_formation in call_data['preconditions']:
            transition_scores.append(1.0)  # Valid transition
        else:
            transition_scores.append(0.0)  # Invalid breaks flow
        
        current_formation = call_data['postconditions'][0]
        formation_trace.append(current_formation)
    
    # Calculate smoothness
    smoothness_score = sum(transition_scores) / len(transition_scores) if transition_scores else 0.0
    
    # Analyze symmetry preservation
    symmetry_trace = [_infer_symmetry(f) for f in formation_trace]
    symmetry_changes = len([i for i in range(1, len(symmetry_trace)) 
                           if symmetry_trace[i] != symmetry_trace[i-1]])
    
    symmetry_analysis = {
        "initial_symmetry": symmetry_trace[0],
        "final_symmetry": symmetry_trace[-1],
        "change_count": symmetry_changes,
        "narrative": "maintains" if symmetry_changes == 0 else
                    "breaks_and_restores" if symmetry_trace[0] == symmetry_trace[-1] else
                    "transforms"
    }
    
    # Spatial coherence
    spatial_coherence = _calculate_spatial_coherence(formation_trace)
    
    # Temporal coherence
    total_beats = sum(CALL_LIBRARY[c]['timing'] for c in call_sequence)
    phrase_aligned = (total_beats % 8 == 0)  # Square dance phrases in 8-beat units
    
    temporal_analysis = {
        "total_beats": total_beats,
        "phrase_aligned": phrase_aligned,
        "sequence_length": len(call_sequence),
        "average_call_timing": total_beats / len(call_sequence)
    }
    
    # Overall flow classification
    if smoothness_score >= 0.9 and spatial_coherence >= 0.8:
        flow_class = "virtuosic_effortless"
    elif smoothness_score >= 0.8:
        flow_class = "seamless_elegant"
    elif smoothness_score >= 0.6:
        flow_class = "smooth_flowing"
    elif smoothness_score >= 0.4:
        flow_class = "adequate_functional"
    else:
        flow_class = "choppy_disjointed"
    
    return {
        "smoothness_score": round(smoothness_score, 2),
        "symmetry_preservation": symmetry_analysis,
        "spatial_coherence": round(spatial_coherence, 2),
        "temporal_coherence": temporal_analysis,
        "flow_classification": flow_class,
        "flow_quality_token": flow_class,  # For compatibility with dimensions
        "interpretation": _interpret_flow_quality(flow_class)
    }

@mcp.tool()
def compare_paths_to_formation(
    target_formation: str,
    path_a: List[str],
    path_b: List[str],
    starting_formation: str = "square_set"
) -> Dict:
    """
    Compare two different call sequences that reach the same formation.
    
    Demonstrates PATH DEPENDENCY and "choreographic hysteresis":
    - Different paths → different dancer momentum
    - Different paths → different symmetry history
    - Different paths → different flow characteristics
    
    Even though final formation is identical!
    """
    # Analyze both paths
    analysis_a = detect_generative_effects(path_a, starting_formation)
    analysis_b = detect_generative_effects(path_b, starting_formation)
    
    flow_a = analyze_flow_quality(path_a, starting_formation)
    flow_b = analyze_flow_quality(path_b, starting_formation)
    
    # Check if both reach target
    final_a = analysis_a['composed_properties']['formation_trace'][-1]
    final_b = analysis_b['composed_properties']['formation_trace'][-1]
    
    both_reach_target = (final_a == target_formation and final_b == target_formation)
    
    # Compare characteristics
    differences = []
    
    # Flow quality difference
    if flow_a['flow_classification'] != flow_b['flow_classification']:
        differences.append({
            "dimension": "flow_quality",
            "path_a": flow_a['flow_classification'],
            "path_b": flow_b['flow_classification'],
            "significance": "Different paths produce different subjective flow experience"
        })
    
    # Symmetry narrative difference
    sym_a = analysis_a['composed_properties']['symmetry_journey']
    sym_b = analysis_b['composed_properties']['symmetry_journey']
    if sym_a != sym_b:
        differences.append({
            "dimension": "symmetry_narrative",
            "path_a": sym_a,
            "path_b": sym_b,
            "significance": "Different symmetry-breaking histories despite same endpoint"
        })
    
    # Timing difference
    timing_a = analysis_a['composed_properties']['total_timing']
    timing_b = analysis_b['composed_properties']['total_timing']
    if timing_a != timing_b:
        differences.append({
            "dimension": "temporal_duration",
            "path_a": timing_a,
            "path_b": timing_b,
            "significance": "Different paths take different time to reach same formation"
        })
    
    # Generative score difference
    gen_a = analysis_a['generative_score']
    gen_b = analysis_b['generative_score']
    differences.append({
        "dimension": "generative_score",
        "path_a": gen_a,
        "path_b": gen_b,
        "difference": abs(gen_a - gen_b),
        "significance": "One path produces more emergent properties than the other"
    })
    
    return {
        "target_formation": target_formation,
        "both_paths_reach_target": both_reach_target,
        "final_formation_a": final_a,
        "final_formation_b": final_b,
        "path_dependency_demonstrated": len(differences) > 0,
        "differences": differences,
        "path_a_analysis": {
            "calls": path_a,
            "generative_score": gen_a,
            "flow_quality": flow_a['flow_classification'],
            "symmetry_narrative": sym_a
        },
        "path_b_analysis": {
            "calls": path_b,
            "generative_score": gen_b,
            "flow_quality": flow_b['flow_classification'],
            "symmetry_narrative": sym_b
        },
        "key_insight": "Same endpoint, different journey - path dependency in discrete morphological systems"
    }

@mcp.tool()
def trace_symmetry_evolution(
    call_sequence: List[str],
    starting_formation: str = "square_set"
) -> Dict:
    """
    Track how symmetry class changes through sequence.
    
    Maps to the "forgetting characteristics" emergence mechanism:
    - Symmetry breaking creates generative opportunities
    - Symmetry restoration is emergent achievement
    - The journey IS the compositional property
    """
    # Validate
    for call in call_sequence:
        if call not in CALL_LIBRARY:
            return {"error": f"Unknown call: {call}"}
    
    # Trace formations
    current_formation = starting_formation
    formation_trace = [starting_formation]
    
    for call in call_sequence:
        call_data = CALL_LIBRARY[call]
        current_formation = call_data['postconditions'][0]
        formation_trace.append(current_formation)
    
    # Trace symmetry states
    symmetry_trace = [_infer_symmetry(f) for f in formation_trace]
    
    # Analyze symmetry narrative
    phases = []
    for i in range(len(symmetry_trace)):
        if i == 0:
            phase = "initial_state"
        elif symmetry_trace[i] == symmetry_trace[i-1]:
            phase = "stable"
        elif symmetry_trace[i] != symmetry_trace[0] and symmetry_trace[i-1] == symmetry_trace[0]:
            phase = "breaking"
        elif symmetry_trace[i] == symmetry_trace[0] and symmetry_trace[i-1] != symmetry_trace[0]:
            phase = "restoring"
        else:
            phase = "transformed"
        
        phases.append({
            "step": i,
            "call": call_sequence[i-1] if i > 0 else "START",
            "formation": formation_trace[i],
            "symmetry": symmetry_trace[i],
            "phase": phase
        })
    
    # Classify overall narrative
    if len(set(symmetry_trace)) == 1:
        narrative_type = "symmetry_preserved"
        generative_potential = "low"
    elif symmetry_trace[0] == symmetry_trace[-1] and len(set(symmetry_trace)) > 1:
        narrative_type = "symmetry_cyclic"
        generative_potential = "high"
    else:
        narrative_type = "symmetry_transformed"
        generative_potential = "medium"
    
    return {
        "call_sequence": call_sequence,
        "symmetry_evolution": phases,
        "symmetry_trace": symmetry_trace,
        "narrative_type": narrative_type,
        "generative_potential": generative_potential,
        "key_transitions": [p for p in phases if p['phase'] in ['breaking', 'restoring']],
        "interpretation": _interpret_symmetry_narrative(narrative_type, phases)
    }

# ============================================================================
# PHASE 2: HIP HOP CASCADE INTEGRATION
# ============================================================================

@mcp.tool()
def translate_hiphop_to_squaredance_params(hiphop_params: Dict[str, str]) -> Dict:
    """
    Translate hip hop rhythmic parameters to square dance choreography.
    
    This is the CROSS-DOMAIN FUNCTOR demonstrating that rhythmic structure
    transcends specific movement vocabularies.
    
    Hip hop rhythmic dimensions map to square dance dimensions via semantic anchors.
    """
    # Extract hip hop parameters (from hip-hop-rhythmic-mcp)
    flow_pattern = hiphop_params.get('flow_pattern', 'moderate_continuous')
    rhythmic_density = hiphop_params.get('rhythmic_density', 'syllables_moderate')
    beat_emphasis = hiphop_params.get('beat_emphasis', 'on_beat')
    breath_pattern = hiphop_params.get('breath_pattern', 'steady_measured')
    delivery_style = hiphop_params.get('delivery_style', 'melodic_sung')
    
    # Map to square dance parameters via semantic anchors
    mapping = {}
    
    # Flow pattern → step pattern + call style
    flow_map = {
        'staccato_choppy': ('quick_allemande', 'syncopated_hash'),
        'moderate_continuous': ('shuffle_step', 'moderate_pace'),
        'legato_flowing': ('walking_step', 'flowing_continuous'),
        'variable_dynamic': ('complex_figure', 'accelerating_figures')
    }
    step, call_style = flow_map.get(flow_pattern, ('shuffle_step', 'moderate_pace'))
    mapping['step_pattern'] = step
    mapping['call_style'] = call_style
    
    # Rhythmic density → timing structure + complexity
    density_map = {
        'syllables_sparse': ('sixteen_beat_figure', 'basic_mainstream'),
        'syllables_moderate': ('thirty_two_beat_chorus', 'plus_level'),
        'syllables_dense': ('sixty_four_beat_sequence', 'advanced_challenging'),
        'syllables_ultra_dense': ('one_twenty_eight_beat_extended', 'c1_complex')
    }
    timing, complexity = density_map.get(rhythmic_density, ('thirty_two_beat_chorus', 'plus_level'))
    mapping['timing_structure'] = timing
    mapping['call_complexity'] = complexity
    
    # Beat emphasis → couple dynamics
    emphasis_map = {
        'on_beat': 'all_eight_circulate',
        'off_beat': 'corner_progression',
        'syncopated': 'heads_sides_alternating',
        'double_time': 'visiting_couples'
    }
    mapping['couple_dynamics'] = emphasis_map.get(beat_emphasis, 'all_eight_circulate')
    
    # Breath pattern → caller interaction
    breath_map = {
        'short_burst': 'competitive_hash',
        'steady_measured': 'standard_calls',
        'extended_phrase': 'patter_calling',
        'variable_adaptive': 'cooperative_flow'
    }
    mapping['caller_interaction'] = breath_map.get(breath_pattern, 'standard_calls')
    
    # Delivery style → formation type
    delivery_map = {
        'aggressive_forceful': 'scatter_promenade',
        'melodic_sung': 'square_set',
        'conversational_relaxed': 'circle_formation',
        'rapid_fire': 'lines_of_four'
    }
    mapping['formation_type'] = delivery_map.get(delivery_style, 'square_set')
    
    # Add defaults for remaining dimensions
    mapping['geometric_symmetry'] = 'four_fold_square'
    mapping['floor_participation'] = 'single_square'
    
    return {
        "hiphop_input": hiphop_params,
        "squaredance_output": mapping,
        "functor_type": "rhythmic_structure_functor",
        "preservation_statement": "Rhythmic structure preserved across movement vocabularies",
        "semantic_anchors_used": [
            "flow_pattern → step_pattern (temporal continuity)",
            "rhythmic_density → timing_structure (complexity scaling)",
            "beat_emphasis → couple_dynamics (accent placement)",
            "breath_pattern → caller_interaction (phrasing logic)",
            "delivery_style → formation_type (energy distribution)"
        ]
    }

# ============================================================================
# PHASE 3: COMPARATIVE ISOMORPHISM VISUALIZATION
# ============================================================================

@mcp.tool()
def compare_choreographic_and_aesthetic_composition(
    dance_sequence: List[str],
    aesthetic_domains: Optional[List[str]] = None,
    domain_results: Optional[List[Dict]] = None,
    starting_formation: str = "square_set"
) -> Dict:
    """
    Demonstrate categorical isomorphism between square dance choreography
    and aesthetic MCP composition.
    
    Two modes:
        1. Simulated: Pass aesthetic_domains as list of domain names (legacy)
        2. Real: Pass domain_results as list of actual MCP server outputs
           Each domain_result should have: domain_id, parameters, vocabulary,
           and optionally optical_properties and metadata.
    
    The client (Claude) gathers results from real MCP servers and passes
    them in. This keeps the server deterministic (Layer 2).
    """
    # Analyze dance sequence
    dance_analysis = detect_generative_effects(dance_sequence, starting_formation)
    dance_flow = analyze_flow_quality(dance_sequence, starting_formation)
    
    # Compose aesthetic domains - real or simulated
    if domain_results and len(domain_results) >= 2:
        aesthetic_composition = _compose_real_domain_results(domain_results)
        aesthetic_domain_names = [dr.get("domain_id", f"domain_{i}") for i, dr in enumerate(domain_results)]
    elif aesthetic_domains:
        aesthetic_composition = _simulate_aesthetic_composition(aesthetic_domains)
        aesthetic_domain_names = aesthetic_domains
    else:
        return {"error": "Provide either domain_results (real) or aesthetic_domains (simulated)"}
    
    is_real = aesthetic_composition.get("source") == "real"
    
    # Extract categorical properties from both
    dance_categorical = {
        "objects": dance_analysis['composed_properties']['formation_trace'],
        "morphisms": dance_sequence,
        "composition_property": dance_flow['flow_classification'],
        "emergence_score": dance_analysis['generative_score'],
        "path_dependency": "demonstrated" if len(dance_sequence) > 1 else "not_applicable",
        "symmetry_breaking": len(set(dance_analysis['composed_properties']['symmetry_narrative'])) > 1
    }
    
    aesthetic_categorical = {
        "objects": aesthetic_composition['domain_states'],
        "morphisms": aesthetic_composition['functors'],
        "composition_property": aesthetic_composition['convergent_features'],
        "emergence_score": aesthetic_composition['generative_score'],
        "path_dependency": "demonstrated",
        "symmetry_breaking": aesthetic_composition['n_conflicts'] > 0 if is_real else aesthetic_composition.get('parameter_conflicts', 0) > 0
    }
    
    # Identify isomorphisms
    isomorphisms = []
    
    # 1. Object-Morphism structure
    n_aesthetic_objects = len(aesthetic_categorical['objects'])
    n_aesthetic_morphisms = len(aesthetic_categorical['morphisms'])
    isomorphisms.append({
        "categorical_concept": "objects_and_morphisms",
        "dance_manifestation": f"{len(dance_categorical['objects'])} formations connected by {len(dance_categorical['morphisms'])} calls",
        "aesthetic_manifestation": f"{n_aesthetic_objects} domain states connected by {n_aesthetic_morphisms} functors",
        "isomorphism": "Both are categories with composable morphisms"
    })
    
    # 2. Compositional emergence
    isomorphisms.append({
        "categorical_concept": "emergent_composition_property",
        "dance_manifestation": f"Flow quality ({dance_categorical['composition_property']}) emerges from call composition",
        "aesthetic_manifestation": f"Convergent features emerge from domain composition",
        "isomorphism": "Properties exist only at compositional level, not in individual morphisms"
    })
    
    # 3. Generative effects
    isomorphisms.append({
        "categorical_concept": "generative_score",
        "dance_manifestation": f"Score: {dance_categorical['emergence_score']}",
        "aesthetic_manifestation": f"Score: {aesthetic_categorical['emergence_score']}",
        "isomorphism": "Both exhibit Φ(s1⊕s2) ≠ Φ(s1)∨Φ(s2) - composition creates new properties"
    })
    
    # 4. Path dependency
    isomorphisms.append({
        "categorical_concept": "path_dependency",
        "dance_manifestation": "Different call sequences to same formation have different flow",
        "aesthetic_manifestation": "Different domain orderings produce different convergent features",
        "isomorphism": "History matters - morphism ORDER affects outcome"
    })
    
    # 5. Symmetry breaking / constraint relaxation
    isomorphisms.append({
        "categorical_concept": "constraint_relaxation",
        "dance_manifestation": "Breaking four-fold symmetry enables new choreographic possibilities",
        "aesthetic_manifestation": "Parameter conflicts enable novel aesthetic blends",
        "isomorphism": "Controlled constraint violation creates generative opportunities"
    })
    
    # Additional isomorphism for real data: structural coupling
    if is_real:
        coupling_types = [f.get("structural_coupling") for f in aesthetic_composition['functors']]
        isomorphisms.append({
            "categorical_concept": "morphism_strength_variation",
            "dance_manifestation": "Some call transitions are smooth, others break flow",
            "aesthetic_manifestation": f"Functor coupling varies: {coupling_types}",
            "isomorphism": "Morphism quality varies - not all compositions are equal",
            "real_data": True
        })
    
    return {
        "dance_system": {
            "sequence": dance_sequence,
            "categorical_properties": dance_categorical,
            "analysis": dance_analysis
        },
        "aesthetic_system": {
            "domains": aesthetic_domain_names,
            "categorical_properties": aesthetic_categorical,
            "composition": aesthetic_composition
        },
        "isomorphisms": isomorphisms,
        "data_source": "real_mcp_servers" if is_real else "simulated",
        "key_theorem": "Square dance choreography and aesthetic MCP composition are categorically isomorphic discrete morphological systems",
        "implications": [
            "Physical grounding validates abstract framework",
            "Decades of square dance practice confirms emergence is real",
            "Same categorical principles govern both domains",
            "Framework generalizes to any compositional system"
        ],
        "visualization_ready": True
    }

# ============================================================================
# PHASE 4: CALL LIBRARY EXPLORER
# ============================================================================

@mcp.tool()
def list_calls(
    level: Optional[str] = None,
    formation: Optional[str] = None,
    movement_type: Optional[str] = None
) -> Dict:
    """
    List available calls with optional filtering.
    
    Args:
        level: Filter by difficulty (basic, plus, advanced)
        formation: Filter by calls usable FROM this formation
        movement_type: Filter by movement type (walk, turn, weave, circle, rotation)
    
    Returns call library with graph connectivity info.
    """
    results = {}
    
    for call_name, call_data in CALL_LIBRARY.items():
        # Apply filters
        if level and call_data.get('level', 'basic') != level:
            continue
        if formation and formation not in call_data['preconditions']:
            continue
        if movement_type and call_data['movement_type'] != movement_type:
            continue
        
        results[call_name] = {
            "timing": call_data['timing'],
            "movement_type": call_data['movement_type'],
            "level": call_data.get('level', 'basic'),
            "preconditions": call_data['preconditions'],
            "postconditions": call_data['postconditions'],
            "symmetry_effect": call_data['symmetry_effect'],
            "geometric_primitive": call_data['geometric_primitive']
        }
    
    # Summary stats
    levels = {}
    for call_data in CALL_LIBRARY.values():
        lv = call_data.get('level', 'basic')
        levels[lv] = levels.get(lv, 0) + 1
    
    return {
        "calls": results,
        "count": len(results),
        "total_in_library": len(CALL_LIBRARY),
        "level_distribution": levels,
        "filters_applied": {
            "level": level,
            "formation": formation,
            "movement_type": movement_type
        }
    }


@mcp.tool()
def get_formation_graph() -> Dict:
    """
    Analyze the formation transition graph defined by the call library.
    
    Returns:
        - nodes: All formations reachable
        - edges: Call-mediated transitions between formations
        - connectivity: Which formations can reach which
        - dead_ends: Formations with no outgoing calls
        - hubs: Formations with most incoming/outgoing connections
    
    This is the categorical structure made explicit: formations are objects,
    calls are morphisms, and this function reveals the category's shape.
    """
    # Build adjacency from call library
    edges = []  # (from_formation, to_formation, call_name, timing)
    formation_outgoing = {}  # formation -> list of (target, call)
    formation_incoming = {}  # formation -> list of (source, call)
    all_formations = set()
    
    for call_name, call_data in CALL_LIBRARY.items():
        for pre in call_data['preconditions']:
            for post in call_data['postconditions']:
                edges.append({
                    "from": pre,
                    "to": post,
                    "call": call_name,
                    "timing": call_data['timing'],
                    "level": call_data.get('level', 'basic'),
                    "symmetry_effect": call_data['symmetry_effect']
                })
                
                all_formations.add(pre)
                all_formations.add(post)
                
                if pre not in formation_outgoing:
                    formation_outgoing[pre] = []
                formation_outgoing[pre].append({"to": post, "call": call_name})
                
                if post not in formation_incoming:
                    formation_incoming[post] = []
                formation_incoming[post].append({"from": pre, "call": call_name})
    
    # Identify dead ends and hubs
    dead_ends = [f for f in all_formations if f not in formation_outgoing]
    
    hub_scores = {}
    for f in all_formations:
        out_count = len(formation_outgoing.get(f, []))
        in_count = len(formation_incoming.get(f, []))
        hub_scores[f] = {"outgoing": out_count, "incoming": in_count, "total": out_count + in_count}
    
    sorted_hubs = sorted(hub_scores.items(), key=lambda x: x[1]['total'], reverse=True)
    
    # Simple reachability (BFS from each formation)
    reachability = {}
    for start in all_formations:
        visited = set()
        queue = [start]
        while queue:
            current = queue.pop(0)
            if current in visited:
                continue
            visited.add(current)
            for edge in formation_outgoing.get(current, []):
                if edge['to'] not in visited:
                    queue.append(edge['to'])
        reachability[start] = sorted(list(visited - {start}))
    
    # Strongly connected components (formations that can reach each other)
    strongly_connected = []
    for f1 in all_formations:
        for f2 in all_formations:
            if f1 < f2 and f2 in reachability.get(f1, []) and f1 in reachability.get(f2, []):
                strongly_connected.append((f1, f2))
    
    return {
        "formations": sorted(list(all_formations)),
        "n_formations": len(all_formations),
        "n_edges": len(edges),
        "edges": edges,
        "hubs": sorted_hubs[:5],
        "dead_ends": dead_ends,
        "reachability": reachability,
        "bidirectional_pairs": strongly_connected,
        "graph_density": round(len(edges) / (len(all_formations) ** 2), 3) if all_formations else 0,
        "categorical_interpretation": {
            "objects": f"{len(all_formations)} formation types",
            "morphisms": f"{len(edges)} call transitions",
            "composable_paths": "Use detect_generative_effects to trace multi-call compositions"
        }
    }


@mcp.tool()
def find_valid_sequences(
    starting_formation: str = "square_set",
    target_formation: Optional[str] = None,
    max_length: int = 4,
    max_results: int = 10,
    level_filter: Optional[str] = None
) -> Dict:
    """
    Find valid call sequences between formations using BFS.
    
    Respects precondition/postcondition constraints. Useful for:
    - Finding all paths between two formations
    - Discovering what sequences are reachable from a starting point
    - Comparing path properties for path dependency analysis
    
    Args:
        starting_formation: Where dancers begin
        target_formation: Where to end up (None = any valid endpoint)
        max_length: Maximum number of calls in sequence
        max_results: Cap on returned sequences
        level_filter: Only use calls at this level or below
    
    Returns valid sequences with basic properties.
    """
    level_order = ['basic', 'plus', 'advanced']
    
    def level_allowed(call_level):
        if not level_filter:
            return True
        try:
            return level_order.index(call_level) <= level_order.index(level_filter)
        except ValueError:
            return True
    
    # BFS for valid sequences
    # State: (current_formation, call_sequence)
    results = []
    queue = [(starting_formation, [])]
    visited_states = set()  # (formation, tuple(sequence)) to avoid loops
    
    while queue and len(results) < max_results:
        current_formation, sequence = queue.pop(0)
        
        if len(sequence) > max_length:
            continue
        
        # Check if this is a valid result
        if len(sequence) > 0:
            if target_formation is None or current_formation == target_formation:
                total_timing = sum(CALL_LIBRARY[c]['timing'] for c in sequence)
                results.append({
                    "sequence": sequence,
                    "final_formation": current_formation,
                    "length": len(sequence),
                    "total_beats": total_timing,
                    "phrase_aligned": total_timing % 8 == 0
                })
                
                if target_formation and current_formation == target_formation:
                    continue  # Found target, don't extend further
        
        if len(sequence) >= max_length:
            continue
        
        # Extend with valid calls
        for call_name, call_data in CALL_LIBRARY.items():
            if not level_allowed(call_data.get('level', 'basic')):
                continue
            if current_formation in call_data['preconditions']:
                new_formation = call_data['postconditions'][0]
                new_sequence = sequence + [call_name]
                state_key = (new_formation, tuple(new_sequence))
                
                if state_key not in visited_states:
                    visited_states.add(state_key)
                    queue.append((new_formation, new_sequence))
    
    # Filter to target if specified
    if target_formation:
        results = [r for r in results if r['final_formation'] == target_formation]
    
    return {
        "starting_formation": starting_formation,
        "target_formation": target_formation or "any",
        "max_length": max_length,
        "level_filter": level_filter or "all",
        "sequences_found": len(results),
        "sequences": results[:max_results],
        "note": f"Capped at {max_results} results" if len(results) >= max_results else "All valid sequences shown"
    }


# ============================================================================
# PHASE 5: STANDALONE AESTHETIC COMPOSITION CONNECTOR
# ============================================================================

@mcp.tool()
def compose_aesthetic_domains(domain_results: List[Dict]) -> Dict:
    """
    Compose real aesthetic domain results into categorical structure.
    
    The client (Claude) gathers outputs from multiple MCP servers
    and passes them here for deterministic categorical analysis.
    
    Each domain_result dict should contain:
        - domain_id: str (e.g. "catastrophe_morph", "diatom_morph", "surface_design")
        - parameters: Dict[str, float] (parameter coordinates from that server)
        - vocabulary: List[str] (visual keywords/tokens from that server)
        - optical_properties: Dict (optional, from servers that provide this)
        - metadata: Dict (optional, any domain-specific data)
    
    Returns:
        Categorical composition analysis including functors between domains,
        convergent features, parameter conflicts, and generative score.
    
    This is the REAL version of aesthetic composition - no simulation.
    Layer 2 deterministic operation, 0 tokens.
    
    Example input (what Claude would gather from two real MCP servers):
        [
            {
                "domain_id": "catastrophe_morph",
                "parameters": {"control_complexity": 0.5, "geometric_sharpness": 0.8,
                               "surface_tension": 0.6, "optical_intensity": 0.7,
                               "aesthetic_intensity": 0.65},
                "vocabulary": ["cusp points", "sharp vertices", "faceted planes"],
                "optical_properties": {"finish": "specular", "refraction": "high"}
            },
            {
                "domain_id": "diatom_morph",
                "parameters": {"radial_symmetry": 0.9, "surface_detail": 0.7,
                               "optical_intensity": 0.5, "pattern_regularity": 0.85},
                "vocabulary": ["radial patterns", "silica frustule", "areolae texture"],
                "optical_properties": {"finish": "iridescent", "transparency": "translucent"}
            }
        ]
    """
    composition = _compose_real_domain_results(domain_results)
    
    if "error" in composition:
        return composition
    
    # Add interpretation layer
    interpretation = {
        "n_domains": composition['n_domains'],
        "coupling_summary": [
            f"{f['from_domain']} → {f['to_domain']}: {f['structural_coupling']} (Jaccard: {f['vocabulary_jaccard']})"
            for f in composition['functors']
        ],
        "emergence_level": (
            "high" if composition['generative_score'] >= 0.6 else
            "moderate" if composition['generative_score'] >= 0.3 else
            "low"
        ),
        "conflict_count": composition['n_conflicts'],
        "convergent_feature_count": len(composition['convergent_features'])
    }
    
    composition['interpretation'] = interpretation
    return composition


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def _infer_symmetry(formation: str) -> str:
    """Infer symmetry class from formation name."""
    if 'square' in formation:
        return 'four_fold'
    elif 'circle' in formation:
        return 'radial'
    elif 'line' in formation:
        return 'bilateral'
    elif 'wave' in formation:
        return 'bilateral'
    elif 'star' in formation:
        return 'radial'
    else:
        return 'asymmetric'

def _apply_symmetry_effect(current: str, effect: str) -> str:
    """Apply symmetry transformation effect."""
    if effect == 'maintains':
        return current
    elif effect == 'breaks':
        return 'asymmetric' if current != 'asymmetric' else current
    elif effect == 'restores':
        return 'four_fold'
    elif effect == 'transforms':
        if current == 'four_fold':
            return 'radial'
        elif current == 'radial':
            return 'bilateral'
        else:
            return current
    return current

def _describe_symmetry_journey(trace: List[str]) -> str:
    """Describe the symmetry narrative arc."""
    if len(set(trace)) == 1:
        return f"maintains_{trace[0]}_throughout"
    elif trace[0] == trace[-1]:
        return f"breaks_{trace[0]}_and_restores"
    else:
        return f"{trace[0]}_to_{trace[-1]}_transformation"

def _calculate_flow_quality(calls: List[str], formations: List[str]) -> str:
    """Calculate flow quality from call-formation sequence."""
    # Check formation compatibility
    valid_transitions = 0
    for i, call in enumerate(calls):
        call_data = CALL_LIBRARY[call]
        if formations[i] in call_data['preconditions']:
            valid_transitions += 1
    
    ratio = valid_transitions / len(calls) if calls else 0
    
    if ratio >= 0.9:
        return "seamless_elegant"
    elif ratio >= 0.7:
        return "smooth_flowing"
    elif ratio >= 0.5:
        return "adequate_functional"
    else:
        return "choppy_disjointed"

def _calculate_spatial_coherence(formations: List[str]) -> float:
    """Calculate spatial coherence from formation trace."""
    # Simple heuristic: fewer formation changes = higher coherence
    unique_formations = len(set(formations))
    total_formations = len(formations)
    
    if total_formations == 1:
        return 1.0
    
    coherence = 1.0 - (unique_formations - 1) / total_formations
    return max(0.0, min(1.0, coherence))

def _interpret_generative_score(score: float) -> str:
    """Interpret generative score value."""
    if score >= 0.7:
        return "High emergence - composition creates substantial new properties"
    elif score >= 0.4:
        return "Moderate emergence - composition shows interesting interactions"
    elif score >= 0.2:
        return "Low emergence - composition is mostly additive"
    else:
        return "Minimal emergence - individual calls dominate"

def _interpret_flow_quality(flow_class: str) -> str:
    """Interpret flow quality classification."""
    interpretations = {
        "virtuosic_effortless": "Exceptional - masterful choreographic composition",
        "seamless_elegant": "Excellent - smooth transitions throughout",
        "smooth_flowing": "Good - generally well-composed sequence",
        "adequate_functional": "Acceptable - some rough transitions",
        "choppy_disjointed": "Poor - multiple formation mismatches"
    }
    return interpretations.get(flow_class, "Unknown flow quality")

def _interpret_symmetry_narrative(narrative_type: str, phases: List[Dict]) -> str:
    """Interpret symmetry evolution narrative."""
    if narrative_type == "symmetry_preserved":
        return "Sequence maintains symmetry - low generative potential but high elegance"
    elif narrative_type == "symmetry_cyclic":
        return "Symmetry breaks and restores - HIGH generative potential, master-level composition"
    else:
        return "Symmetry transforms - medium generative potential, exploratory choreography"

def _simulate_aesthetic_composition(domains: List[str]) -> Dict:
    """
    LEGACY fallback: Simulate aesthetic domain composition.
    Used only when real domain results are not provided.
    """
    return {
        "domain_states": [f"{d}_parameters" for d in domains],
        "functors": [f"{domains[i]}_to_{domains[i+1]}_functor" for i in range(len(domains)-1)],
        "convergent_features": "color_temperature_convergence",
        "generative_score": 0.65,
        "parameter_conflicts": 2,
        "resolution_strategy": "hierarchical_authority",
        "source": "simulated"
    }


def _compose_real_domain_results(domain_results: List[Dict]) -> Dict:
    """
    Compose REAL aesthetic domain results into categorical structure.
    
    Each domain_result should have:
        - domain_id: str (e.g. "catastrophe_morph", "diatom_morph")
        - parameters: Dict[str, float] (5D coordinates or equivalent)
        - vocabulary: List[str] (visual keywords)
        - optical_properties: Dict (optional)
        - metadata: Dict (optional - anything domain-specific)
    
    Returns categorical composition analysis (deterministic, 0 tokens).
    """
    if len(domain_results) < 2:
        return {
            "error": "Need at least 2 domains for composition",
            "source": "real"
        }
    
    # Extract domain states (objects in the category)
    domain_states = []
    for dr in domain_results:
        domain_states.append({
            "domain_id": dr.get("domain_id", "unknown"),
            "parameter_count": len(dr.get("parameters", {})),
            "vocabulary_size": len(dr.get("vocabulary", [])),
            "has_optical": "optical_properties" in dr and dr["optical_properties"] is not None
        })
    
    # Build functors (morphisms between adjacent domains)
    functors = []
    for i in range(len(domain_results) - 1):
        d1 = domain_results[i]
        d2 = domain_results[i + 1]
        
        p1 = d1.get("parameters", {})
        p2 = d2.get("parameters", {})
        
        # Find shared parameter names (structural overlap)
        shared_params = set(p1.keys()) & set(p2.keys())
        
        # Calculate parameter distance for shared dimensions
        if shared_params:
            distances = {k: abs(p1[k] - p2[k]) for k in shared_params}
            avg_distance = sum(distances.values()) / len(distances)
        else:
            distances = {}
            avg_distance = None
        
        # Vocabulary overlap (convergent features)
        v1 = set(d1.get("vocabulary", []))
        v2 = set(d2.get("vocabulary", []))
        vocab_overlap = v1 & v2
        vocab_union = v1 | v2
        jaccard = len(vocab_overlap) / len(vocab_union) if vocab_union else 0.0
        
        functors.append({
            "from_domain": d1.get("domain_id", f"domain_{i}"),
            "to_domain": d2.get("domain_id", f"domain_{i+1}"),
            "shared_parameters": list(shared_params),
            "parameter_distances": distances,
            "avg_parameter_distance": round(avg_distance, 4) if avg_distance is not None else None,
            "vocabulary_overlap": list(vocab_overlap)[:10],  # Cap for readability
            "vocabulary_jaccard": round(jaccard, 4),
            "structural_coupling": "tight" if jaccard > 0.3 else "moderate" if jaccard > 0.1 else "loose"
        })
    
    # Detect convergent features (properties emerging from composition)
    all_vocabs = [set(dr.get("vocabulary", [])) for dr in domain_results]
    if len(all_vocabs) >= 2:
        # Words appearing in 2+ domains
        from collections import Counter
        word_counts = Counter()
        for v in all_vocabs:
            word_counts.update(v)
        convergent_features = [w for w, c in word_counts.items() if c >= 2]
    else:
        convergent_features = []
    
    # Detect parameter conflicts
    all_params = {}
    for dr in domain_results:
        for k, v in dr.get("parameters", {}).items():
            if k not in all_params:
                all_params[k] = []
            all_params[k].append((dr.get("domain_id", "?"), v))
    
    conflicts = []
    for param, values in all_params.items():
        if len(values) >= 2:
            vals = [v for _, v in values]
            spread = max(vals) - min(vals)
            if spread > 0.3:  # Threshold for meaningful conflict
                conflicts.append({
                    "parameter": param,
                    "spread": round(spread, 4),
                    "domain_values": {did: round(v, 4) for did, v in values}
                })
    
    # Calculate generative score
    # Based on: functor diversity + convergent features + productive conflicts
    functor_diversity = len(set(f["structural_coupling"] for f in functors)) / 3.0
    convergent_ratio = min(len(convergent_features) / 10.0, 1.0)
    conflict_ratio = min(len(conflicts) / 5.0, 1.0)
    
    generative_score = round(
        0.3 * functor_diversity + 
        0.4 * convergent_ratio + 
        0.3 * conflict_ratio,
        3
    )
    
    return {
        "domain_states": domain_states,
        "functors": functors,
        "convergent_features": convergent_features[:20],
        "parameter_conflicts": conflicts,
        "n_conflicts": len(conflicts),
        "generative_score": generative_score,
        "resolution_strategy": "empirical_composition",
        "source": "real",
        "n_domains": len(domain_results)
    }

# ============================================================================
# PHASE 2.6: RHYTHMIC PRESETS — Temporal Composition
# ============================================================================
#
# Square dancing has a natural 5D aesthetic parameter space derived from
# the choreographic properties that emerge during execution. These are NOT
# the olog dimensions (which are categorical tokens) but continuous
# coordinates capturing the *feel* of a choreographic moment.
#
# This enables the Lushy Aesthetic Dynamics framework to treat square
# dance choreography as a composable aesthetic domain alongside
# microscopy, catastrophe theory, diatom morphology, etc.
# ============================================================================



# --- Normalized Parameter Space [0.0, 1.0] ---

SQUAREDANCE_PARAMETER_NAMES = [
    "kinetic_energy",       # 0.0 = standing still / bow, 1.0 = vigorous swing/hash
    "spatial_complexity",   # 0.0 = simple square set, 1.0 = complex weave/chain
    "rotational_momentum",  # 0.0 = linear pass-thru, 1.0 = full star/swing rotation
    "group_cohesion",       # 0.0 = scattered/individual, 1.0 = locked formation unison
    "rhythmic_regularity",  # 0.0 = syncopated patter hash, 1.0 = metronomic singing call
]

# --- Canonical Choreographic States ---
# Each state represents a recognizable "moment" in square dance with
# distinct aesthetic character, mapped to the 5D parameter space.

SQUAREDANCE_COORDS = {
    "opening_promenade": {
        "kinetic_energy": 0.45,
        "spatial_complexity": 0.20,
        "rotational_momentum": 0.30,
        "group_cohesion": 0.85,
        "rhythmic_regularity": 0.90,
    },
    "hash_breakdown": {
        "kinetic_energy": 0.90,
        "spatial_complexity": 0.85,
        "rotational_momentum": 0.60,
        "group_cohesion": 0.40,
        "rhythmic_regularity": 0.15,
    },
    "singing_call_figure": {
        "kinetic_energy": 0.55,
        "spatial_complexity": 0.50,
        "rotational_momentum": 0.45,
        "group_cohesion": 0.75,
        "rhythmic_regularity": 0.95,
    },
    "grand_right_left_weave": {
        "kinetic_energy": 0.70,
        "spatial_complexity": 0.90,
        "rotational_momentum": 0.80,
        "group_cohesion": 0.55,
        "rhythmic_regularity": 0.60,
    },
    "star_turn": {
        "kinetic_energy": 0.60,
        "spatial_complexity": 0.35,
        "rotational_momentum": 0.95,
        "group_cohesion": 0.90,
        "rhythmic_regularity": 0.70,
    },
    "ocean_wave_flow": {
        "kinetic_energy": 0.50,
        "spatial_complexity": 0.65,
        "rotational_momentum": 0.40,
        "group_cohesion": 0.70,
        "rhythmic_regularity": 0.55,
    },
    "scatter_promenade": {
        "kinetic_energy": 0.75,
        "spatial_complexity": 0.45,
        "rotational_momentum": 0.25,
        "group_cohesion": 0.20,
        "rhythmic_regularity": 0.35,
    },
}

# --- Phase 2.6 Rhythmic Presets ---
# Each preset defines a temporal oscillation between two canonical states.
# Periods chosen to create useful resonances with other Lushy domains:
#   microscopy [10, 16, 20, 24, 30], nuclear [15, 18],
#   catastrophe [15, 18, 20, 22, 25], diatom [12, 15, 18, 20, 30]

SQUAREDANCE_RHYTHMIC_PRESETS = {
    "energy_cycle": {
        "state_a": "opening_promenade",
        "state_b": "hash_breakdown",
        "pattern": "sinusoidal",
        "num_cycles": 3,
        "steps_per_cycle": 20,
        "description": "Smooth cycling between calm promenade and energetic hash"
    },
    "formation_sweep": {
        "state_a": "star_turn",
        "state_b": "grand_right_left_weave",
        "pattern": "triangular",
        "num_cycles": 2,
        "steps_per_cycle": 24,
        "description": "Linear ramp between tight star and expansive grand chain"
    },
    "caller_intensity": {
        "state_a": "singing_call_figure",
        "state_b": "hash_breakdown",
        "pattern": "sinusoidal",
        "num_cycles": 4,
        "steps_per_cycle": 16,
        "description": "Pulse between structured singing call and freeform hash"
    },
    "cohesion_pulse": {
        "state_a": "star_turn",
        "state_b": "scatter_promenade",
        "pattern": "sinusoidal",
        "num_cycles": 4,
        "steps_per_cycle": 15,
        "description": "Breathing between tight locked formation and scattered movement"
    },
    "complexity_wave": {
        "state_a": "opening_promenade",
        "state_b": "ocean_wave_flow",
        "pattern": "sinusoidal",
        "num_cycles": 3,
        "steps_per_cycle": 18,
        "description": "Gentle oscillation from simple beginnings to flowing waves"
    },
}


def _generate_squaredance_oscillation(
    num_steps: int,
    num_cycles: float,
    pattern: str,
) -> List[float]:
    """Generate oscillation pattern values in [0, 1]."""
    result = []
    for i in range(num_steps):
        t = 2.0 * math.pi * num_cycles * i / num_steps
        if pattern == "sinusoidal":
            result.append(0.5 * (1.0 + math.sin(t)))
        elif pattern == "triangular":
            t_norm = (t / (2.0 * math.pi)) % 1.0
            result.append(2.0 * t_norm if t_norm < 0.5 else 2.0 * (1.0 - t_norm))
        elif pattern == "square":
            t_norm = (t / (2.0 * math.pi)) % 1.0
            result.append(0.0 if t_norm < 0.5 else 1.0)
        else:
            result.append(0.5)
    return result


def _generate_squaredance_preset_trajectory(preset_config: dict) -> List[Dict[str, float]]:
    """Generate a full preset trajectory as list of state dicts."""
    state_a = SQUAREDANCE_COORDS[preset_config["state_a"]]
    state_b = SQUAREDANCE_COORDS[preset_config["state_b"]]
    total_steps = preset_config["num_cycles"] * preset_config["steps_per_cycle"]
    alphas = _generate_squaredance_oscillation(
        total_steps, preset_config["num_cycles"], preset_config["pattern"]
    )
    trajectory = []
    for alpha in alphas:
        state = {}
        for p in SQUAREDANCE_PARAMETER_NAMES:
            state[p] = round(state_a[p] * (1.0 - alpha) + state_b[p] * alpha, 6)
        trajectory.append(state)
    return trajectory


# --- Phase 2.6 Tools ---

@mcp.tool()
def get_squaredance_rhythmic_presets() -> Dict:
    """
    List all Phase 2.6 rhythmic presets for square dance choreography.

    Returns preset names, periods, patterns, state endpoints, and
    descriptions. These presets define temporal oscillations between
    canonical choreographic states in 5D parameter space.

    Periods: [15, 16, 18, 20, 24] — chosen for resonance with
    microscopy, nuclear, catastrophe, and diatom domain periods.

    Layer 1: Pure taxonomy lookup (0 tokens).
    """
    result = {}
    for name, cfg in SQUAREDANCE_RHYTHMIC_PRESETS.items():
        result[name] = {
            "period": cfg["steps_per_cycle"],
            "total_steps": cfg["num_cycles"] * cfg["steps_per_cycle"],
            "pattern": cfg["pattern"],
            "state_a": cfg["state_a"],
            "state_b": cfg["state_b"],
            "description": cfg["description"],
        }
    return {
        "domain": "squaredance",
        "presets": result,
        "parameter_names": SQUAREDANCE_PARAMETER_NAMES,
        "n_canonical_states": len(SQUAREDANCE_COORDS),
        "periods": sorted(set(cfg["steps_per_cycle"] for cfg in SQUAREDANCE_RHYTHMIC_PRESETS.values())),
    }


@mcp.tool()
def get_squaredance_coordinates() -> Dict:
    """
    Get all canonical choreographic state coordinates in 5D parameter space.

    Each state is a recognizable square dance moment with distinct
    aesthetic character, mapped to normalized [0, 1] coordinates.

    Layer 1: Pure taxonomy lookup (0 tokens).
    """
    return {
        "domain": "squaredance",
        "parameter_names": SQUAREDANCE_PARAMETER_NAMES,
        "states": SQUAREDANCE_COORDS,
        "n_states": len(SQUAREDANCE_COORDS),
        "parameter_semantics": {
            "kinetic_energy": "0.0 = standing still / bow, 1.0 = vigorous swing/hash",
            "spatial_complexity": "0.0 = simple square set, 1.0 = complex weave/chain",
            "rotational_momentum": "0.0 = linear pass-thru, 1.0 = full star/swing rotation",
            "group_cohesion": "0.0 = scattered/individual, 1.0 = locked formation unison",
            "rhythmic_regularity": "0.0 = syncopated patter hash, 1.0 = metronomic singing call",
        },
    }


@mcp.tool()
def generate_squaredance_rhythmic_sequence(
    state_a_id: str,
    state_b_id: str,
    oscillation_pattern: str = "sinusoidal",
    num_cycles: int = 3,
    steps_per_cycle: int = 20,
    phase_offset: float = 0.0,
) -> Dict:
    """
    Generate rhythmic oscillation between two choreographic states.

    Layer 2: Temporal composition (0 tokens).

    Creates periodic transitions cycling between canonical square dance
    states in 5D parameter space. Useful for composing with other
    Lushy domains via forced orbit integration.

    Args:
        state_a_id: Starting state (opening_promenade, hash_breakdown, etc.)
        state_b_id: Alternating state
        oscillation_pattern: "sinusoidal", "triangular", or "square"
        num_cycles: Number of complete A→B→A cycles
        steps_per_cycle: Samples per cycle (= period for limit cycle detection)
        phase_offset: Starting phase (0.0 = A, 0.5 = B)

    Returns:
        Sequence with parameter states at each step, pattern info, and
        phase points for integration with aesthetic-dynamics-core.
    """
    if state_a_id not in SQUAREDANCE_COORDS:
        return {"error": f"Unknown state: {state_a_id}. Valid: {list(SQUAREDANCE_COORDS.keys())}"}
    if state_b_id not in SQUAREDANCE_COORDS:
        return {"error": f"Unknown state: {state_b_id}. Valid: {list(SQUAREDANCE_COORDS.keys())}"}
    if oscillation_pattern not in ("sinusoidal", "triangular", "square"):
        return {"error": f"Unknown pattern: {oscillation_pattern}"}

    state_a = SQUAREDANCE_COORDS[state_a_id]
    state_b = SQUAREDANCE_COORDS[state_b_id]
    total_steps = num_cycles * steps_per_cycle

    alphas = _generate_squaredance_oscillation(total_steps, num_cycles, oscillation_pattern)

    # Apply phase offset
    if phase_offset > 0.0:
        offset_steps = int(phase_offset * steps_per_cycle)
        alphas = alphas[offset_steps:] + alphas[:offset_steps]

    sequence = []
    for step, alpha in enumerate(alphas):
        state = {}
        for p in SQUAREDANCE_PARAMETER_NAMES:
            state[p] = round(state_a[p] * (1.0 - alpha) + state_b[p] * alpha, 6)
        sequence.append({
            "step": step,
            "phase": round(alpha, 6),
            "state": state,
        })

    return {
        "domain": "squaredance",
        "state_a": state_a_id,
        "state_b": state_b_id,
        "pattern": oscillation_pattern,
        "num_cycles": num_cycles,
        "steps_per_cycle": steps_per_cycle,
        "total_steps": total_steps,
        "phase_offset": phase_offset,
        "parameter_names": SQUAREDANCE_PARAMETER_NAMES,
        "sequence": sequence,
    }


@mcp.tool()
def apply_squaredance_rhythmic_preset(preset_name: str) -> Dict:
    """
    Apply a curated Phase 2.6 rhythmic preset.

    Layer 2: Deterministic sequence generation (0 tokens).

    Available presets:
        energy_cycle (20): promenade ↔ hash (energy sweep)
        formation_sweep (24): star ↔ grand chain (spatial ramp)
        caller_intensity (16): singing call ↔ hash (structure pulse)
        cohesion_pulse (15): star ↔ scatter (cohesion breathing)
        complexity_wave (18): promenade ↔ ocean wave (complexity oscillation)

    Returns complete oscillation sequence with parameter states at each step.
    """
    if preset_name not in SQUAREDANCE_RHYTHMIC_PRESETS:
        return {
            "error": f"Unknown preset: {preset_name}",
            "available": list(SQUAREDANCE_RHYTHMIC_PRESETS.keys()),
        }

    cfg = SQUAREDANCE_RHYTHMIC_PRESETS[preset_name]
    trajectory = _generate_squaredance_preset_trajectory(cfg)

    return {
        "domain": "squaredance",
        "preset_name": preset_name,
        "period": cfg["steps_per_cycle"],
        "total_steps": len(trajectory),
        "pattern": cfg["pattern"],
        "state_a": cfg["state_a"],
        "state_b": cfg["state_b"],
        "description": cfg["description"],
        "parameter_names": SQUAREDANCE_PARAMETER_NAMES,
        "sequence": [
            {"step": i, "state": s} for i, s in enumerate(trajectory)
        ],
    }


# ============================================================================
# PHASE 2.7: ATTRACTOR VISUALIZATION — Prompt Generation
# ============================================================================
#
# Maps choreographic parameter states to visual vocabulary suitable for
# image generation (ComfyUI, Stable Diffusion, DALL-E, etc.).
#
# Visual types represent distinct aesthetic renderings of square dance
# scenes. Nearest-neighbor matching in 5D space selects the best
# vocabulary for any parameter coordinate.
# ============================================================================

SQUAREDANCE_VISUAL_TYPES = {
    "barn_dance_warm": {
        "coords": {
            "kinetic_energy": 0.45,
            "spatial_complexity": 0.25,
            "rotational_momentum": 0.35,
            "group_cohesion": 0.85,
            "rhythmic_regularity": 0.85,
        },
        "keywords": [
            "warm barn interior with golden lantern light",
            "couples in unified promenade circle",
            "wooden dance floor with motion blur at feet",
            "communal intimacy and folk tradition",
            "soft amber glow on cotton and denim fabrics",
            "four-fold square symmetry visible from above",
            "relaxed rhythmic sway captured mid-stride",
        ],
    },
    "competition_precise": {
        "coords": {
            "kinetic_energy": 0.80,
            "spatial_complexity": 0.75,
            "rotational_momentum": 0.70,
            "group_cohesion": 0.90,
            "rhythmic_regularity": 0.80,
        },
        "keywords": [
            "crisp geometric formation under bright stage lighting",
            "dancers in matched western attire with sharp creases",
            "precise angular arm positions at exact 90-degree extensions",
            "polished hardwood floor reflecting overhead spotlights",
            "synchronized eight-dancer figure with military precision",
            "clean diagonal sight-lines and radial star geometry",
            "high-contrast lighting with deep cast shadows",
        ],
    },
    "folk_festival_vibrant": {
        "coords": {
            "kinetic_energy": 0.85,
            "spatial_complexity": 0.60,
            "rotational_momentum": 0.55,
            "group_cohesion": 0.50,
            "rhythmic_regularity": 0.40,
        },
        "keywords": [
            "outdoor festival with colored bunting and string lights",
            "dancers in bright patterned skirts mid-twirl",
            "energetic scatter with multiple couples in different phases",
            "dust or confetti particles suspended in evening light",
            "joyful kinetic blur suggesting rapid movement",
            "mixed formations across wide grassy dance area",
            "warm sunset backlighting silhouetting raised arms",
        ],
    },
    "diagram_technical": {
        "coords": {
            "kinetic_energy": 0.30,
            "spatial_complexity": 0.90,
            "rotational_momentum": 0.50,
            "group_cohesion": 0.65,
            "rhythmic_regularity": 0.70,
        },
        "keywords": [
            "clean top-down choreography diagram on cream paper",
            "numbered dancer positions with directional flow arrows",
            "formation geometry drawn with precise compass arcs",
            "call notation labels in serif typography",
            "dotted path traces showing weave and chain routing",
            "bilateral symmetry axis marked in fine red line",
            "instructional clarity with minimal decoration",
        ],
    },
    "twilight_romantic": {
        "coords": {
            "kinetic_energy": 0.40,
            "spatial_complexity": 0.30,
            "rotational_momentum": 0.80,
            "group_cohesion": 0.80,
            "rhythmic_regularity": 0.75,
        },
        "keywords": [
            "couple swing under string lights at dusk",
            "soft bokeh background of distant barn and trees",
            "flowing skirt fabric captured in rotational arc",
            "intimate partner connection with joined hands",
            "cool blue twilight sky gradient above warm ground light",
            "slow-shutter rotational blur around sharp central faces",
            "romantic pastoral atmosphere with firefly points of light",
        ],
    },
}


def _extract_squaredance_visual_vocabulary(
    state: Dict[str, float],
    strength: float = 1.0,
) -> Dict:
    """
    Map a 5D parameter state to nearest visual type and return keywords.

    Pure Layer 2 deterministic computation (0 tokens).
    """
    # Compute distance to each visual type
    best_type = None
    best_dist = float("inf")
    for vtype, vdata in SQUAREDANCE_VISUAL_TYPES.items():
        dist_sq = 0.0
        for p in SQUAREDANCE_PARAMETER_NAMES:
            diff = state.get(p, 0.5) - vdata["coords"][p]
            dist_sq += diff * diff
        dist = math.sqrt(dist_sq)
        if dist < best_dist:
            best_dist = dist
            best_type = vtype

    vdata = SQUAREDANCE_VISUAL_TYPES[best_type]

    # Weight keywords by strength
    if strength >= 0.7:
        keywords = vdata["keywords"]
    elif strength >= 0.4:
        keywords = vdata["keywords"][:5]
    else:
        keywords = vdata["keywords"][:3]

    return {
        "nearest_type": best_type,
        "distance": round(best_dist, 4),
        "keywords": keywords,
        "strength": strength,
        "coords_used": {p: round(state.get(p, 0.5), 4) for p in SQUAREDANCE_PARAMETER_NAMES},
    }


@mcp.tool()
def extract_squaredance_visual_vocabulary(
    state: Optional[Dict[str, float]] = None,
    state_id: Optional[str] = None,
    strength: float = 1.0,
) -> Dict:
    """
    Extract visual vocabulary from square dance parameter coordinates.

    Layer 2: Deterministic vocabulary mapping (0 tokens).

    Maps a 5D parameter state to the nearest canonical visual type
    and returns image-generation-ready keywords. Provide either a
    raw state dict or a canonical state_id.

    Args:
        state: Parameter coordinates dict (kinetic_energy, spatial_complexity, etc.)
        state_id: Canonical state name (opening_promenade, hash_breakdown, etc.)
        strength: Keyword weight multiplier [0.0, 1.0]

    Returns:
        Nearest visual type, keywords, distance, and coordinate echo.
    """
    if state_id and state_id in SQUAREDANCE_COORDS:
        state = SQUAREDANCE_COORDS[state_id]
    elif state is None:
        return {"error": "Provide either 'state' dict or 'state_id'"}

    return _extract_squaredance_visual_vocabulary(state, strength)


@mcp.tool()
def generate_squaredance_visualization_prompt(
    preset_name: str = "",
    state_id: str = "",
    custom_state: Optional[Dict[str, float]] = None,
    mode: str = "composite",
    style_modifier: str = "",
    keyframe_count: int = 4,
) -> Dict:
    """
    Generate image-generation prompt from choreographic state or preset.

    Layer 2: Deterministic prompt synthesis (0 tokens).

    Translates square dance parameter coordinates into visual prompts
    suitable for ComfyUI, Stable Diffusion, DALL-E, etc.

    Modes:
        composite: Single blended prompt from state
        sequence: Multiple keyframe prompts from preset trajectory

    Args:
        preset_name: Phase 2.6 preset (energy_cycle, formation_sweep, etc.)
        state_id: Canonical state for composite mode
        custom_state: Optional custom 5D coordinates
        mode: "composite" or "sequence"
        style_modifier: Optional prefix ("oil painting", "watercolor", etc.)
        keyframe_count: Number of keyframes for sequence mode (default: 4)

    Returns:
        Prompt string(s) with vocabulary details and state metadata.
    """
    if mode == "sequence":
        # Generate keyframe prompts from preset
        if not preset_name or preset_name not in SQUAREDANCE_RHYTHMIC_PRESETS:
            return {
                "error": f"Sequence mode requires valid preset_name. Available: {list(SQUAREDANCE_RHYTHMIC_PRESETS.keys())}",
            }
        cfg = SQUAREDANCE_RHYTHMIC_PRESETS[preset_name]
        trajectory = _generate_squaredance_preset_trajectory(cfg)
        total = len(trajectory)

        # Extract evenly spaced keyframes
        keyframes = []
        for k in range(keyframe_count):
            idx = int(k * total / keyframe_count)
            state = trajectory[idx]
            vocab = _extract_squaredance_visual_vocabulary(state)
            prompt_parts = []
            if style_modifier:
                prompt_parts.append(style_modifier)
            prompt_parts.extend(vocab["keywords"])
            keyframes.append({
                "keyframe": k,
                "step": idx,
                "prompt": ", ".join(prompt_parts),
                "nearest_type": vocab["nearest_type"],
                "distance": vocab["distance"],
                "state": {p: round(state[p], 4) for p in SQUAREDANCE_PARAMETER_NAMES},
            })

        return {
            "mode": "sequence",
            "preset": preset_name,
            "description": cfg["description"],
            "period": cfg["steps_per_cycle"],
            "keyframe_count": keyframe_count,
            "keyframes": keyframes,
        }

    else:  # composite
        # Determine state
        if custom_state:
            state = custom_state
        elif state_id and state_id in SQUAREDANCE_COORDS:
            state = SQUAREDANCE_COORDS[state_id]
        elif preset_name and preset_name in SQUAREDANCE_RHYTHMIC_PRESETS:
            # Use midpoint of preset as representative state
            cfg = SQUAREDANCE_RHYTHMIC_PRESETS[preset_name]
            sa = SQUAREDANCE_COORDS[cfg["state_a"]]
            sb = SQUAREDANCE_COORDS[cfg["state_b"]]
            state = {p: (sa[p] + sb[p]) / 2.0 for p in SQUAREDANCE_PARAMETER_NAMES}
        else:
            return {
                "error": "Composite mode requires state_id, custom_state, or preset_name",
                "available_states": list(SQUAREDANCE_COORDS.keys()),
                "available_presets": list(SQUAREDANCE_RHYTHMIC_PRESETS.keys()),
            }

        vocab = _extract_squaredance_visual_vocabulary(state)
        prompt_parts = []
        if style_modifier:
            prompt_parts.append(style_modifier)
        prompt_parts.extend(vocab["keywords"])

        return {
            "mode": "composite",
            "prompt": ", ".join(prompt_parts),
            "nearest_type": vocab["nearest_type"],
            "distance": vocab["distance"],
            "state": {p: round(state.get(p, 0.5), 4) for p in SQUAREDANCE_PARAMETER_NAMES},
            "vocabulary": vocab,
        }


@mcp.tool()
def get_squaredance_visual_types() -> Dict:
    """
    List all canonical visual types for square dance choreography.

    Layer 1: Pure taxonomy lookup (0 tokens).

    Returns 5 visual types spanning the choreographic morphospace:
        barn_dance_warm — Communal warmth, golden light, folk tradition
        competition_precise — Geometric precision, stage lighting, discipline
        folk_festival_vibrant — Outdoor energy, color, kinetic joy
        diagram_technical — Instructional clarity, formation geometry
        twilight_romantic — Intimate rotation, dusk atmosphere, bokeh
    """
    result = {}
    for vtype, vdata in SQUAREDANCE_VISUAL_TYPES.items():
        result[vtype] = {
            "coords": vdata["coords"],
            "keywords": vdata["keywords"],
            "keyword_count": len(vdata["keywords"]),
        }
    return {
        "domain": "squaredance",
        "visual_types": result,
        "n_types": len(result),
        "parameter_names": SQUAREDANCE_PARAMETER_NAMES,
    }


# ============================================================================
# SERVER INFO UPDATE — Advertise Phase 2.6 + 2.7 Capabilities
# ============================================================================

@mcp.tool()
def get_server_info() -> Dict:
    """
    Get comprehensive information about the Square Dance Choreography MCP server.

    Includes Phase 2.6 rhythmic preset and Phase 2.7 attractor visualization
    capabilities for integration with Lushy Aesthetic Dynamics framework.
    """
    return {
        "server_name": "squaredance-choreography-extended",
        "version": "2.6.0",
        "description": (
            "Square dancing as a categorically composable aesthetic domain. "
            "Provides call sequence analysis, generative effects detection, "
            "cross-domain functor translation, and Phase 2.6/2.7 rhythmic "
            "preset + attractor visualization support for multi-domain "
            "aesthetic composition."
        ),
        "capabilities": {
            "generative_effects_analysis": True,
            "flow_quality_measurement": True,
            "symmetry_evolution_tracking": True,
            "path_dependency_comparison": True,
            "cross_domain_translation": True,
            "formation_graph_analysis": True,
            "aesthetic_domain_composition": True,
        },
        "phase_2_6_enhancements": {
            "rhythmic_presets": True,
            "n_presets": len(SQUAREDANCE_RHYTHMIC_PRESETS),
            "presets": {
                name: {
                    "period": cfg["steps_per_cycle"],
                    "pattern": cfg["pattern"],
                    "states": f"{cfg['state_a']} ↔ {cfg['state_b']}",
                }
                for name, cfg in SQUAREDANCE_RHYTHMIC_PRESETS.items()
            },
            "periods": sorted(set(
                cfg["steps_per_cycle"] for cfg in SQUAREDANCE_RHYTHMIC_PRESETS.values()
            )),
            "parameter_space": {
                "dimensions": 5,
                "parameter_names": SQUAREDANCE_PARAMETER_NAMES,
                "canonical_states": len(SQUAREDANCE_COORDS),
            },
        },
        "phase_2_7_enhancements": {
            "attractor_visualization": True,
            "visual_types": list(SQUAREDANCE_VISUAL_TYPES.keys()),
            "n_visual_types": len(SQUAREDANCE_VISUAL_TYPES),
            "prompt_modes": ["composite", "sequence"],
            "supported_generators": ["ComfyUI", "Stable Diffusion", "DALL-E", "Midjourney"],
        },
        "integration_with_lushy": {
            "domain_id": "squaredance",
            "composable_with": [
                "microscopy", "nuclear", "catastrophe", "diatom",
                "heraldic", "surface_design",
            ],
            "shared_periods": {
                "microscopy": [16, 20, 24],
                "nuclear": [15, 18],
                "catastrophe": [15, 16, 18, 20],
                "diatom": [15, 18, 20],
            },
            "expected_lcm_hubs": [60, 120, 360],
            "novel_attractor_potential": "high — 5 periods with broad coverage",
        },
        "call_library": {
            "total_calls": len(CALL_LIBRARY),
            "levels": ["basic", "plus", "advanced"],
        },
        "layer_costs": {
            "taxonomy_lookups": "0 tokens (Layer 1)",
            "rhythmic_sequences": "0 tokens (Layer 2)",
            "visual_vocabulary": "0 tokens (Layer 2)",
            "prompt_generation": "0 tokens (Layer 2)",
        },
    }


# Run server
if __name__ == "__main__":
    mcp.run()
