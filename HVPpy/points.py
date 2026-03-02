
threshold = {
    1: 16, 2: 16, 3:16,
    4: 4, 5: 4, 6: 4
    }
problematic = {
    1: set(), 2: {16}, 3: {0,4,16}, # 0 and 4 are removable singularities but they're a pain to remove
    4: {0,16}, 5: {0,4}, 6: {0,4}      # Likewise, 0 (and 4 for E5) are removable
    }
problematic_with_theta = {
    1: {0,4,16}, 2: {0,4,16}, 3: {0,4,16},
    4: {0,4,16}, 5: {0,4,16}, 6: {0,4,16}
    }
