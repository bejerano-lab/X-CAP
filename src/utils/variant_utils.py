class Variant():
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom.replace('chr', '')
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt
    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.chrom == other.chrom and \
                   self.pos == other.pos and \
                   self.ref == other.ref and \
                   self.alt == other.alt
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)
   
    def __hash__(self):
        return hash((self.chrom, self.pos, self.ref, self.alt))

    def __str__(self):
        return "Chr: {}, pos: {}, ref: {}, alt: {}".format(
            self.chrom, self.pos, self.ref, self.alt)
