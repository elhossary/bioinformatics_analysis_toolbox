from numpy import diff, where, split


class PatternLocationsFetcher:

    def __init__(self, seq_str, base, orientation, max_interruption):
        self.seq_str = seq_str
        self.base = base
        self.orientation = orientation
        self.max_interruption = max_interruption

    def fetch_locations(self):
        complement_base = lambda x: "T" if self.base == "A" else "A"
        indices = []
        if self.orientation == 'f':
            indices = [i for i, a in enumerate(self.seq_str, 1) if a == self.base]
        if self.orientation == 'r':
            indices = [i for i, a in enumerate(self.seq_str, 1) if a == complement_base(self.base)]
        indices.sort()
        # Get all signals to any length with max interruption
        pattern_locations = list(map(list, split(indices, where(diff(indices) > self.max_interruption)[0] + 1)))
        return pattern_locations

