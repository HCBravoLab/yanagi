
# A Segment Key that is used to identify a segment in the segments graph (SG)
class SG_Key:
    exs = []
    exs_str = ""
    pos = 0

    def __init__(self):
        self.exs = []
        self.exs_str = ""
        self.pos = 0

    def addEx(self, ex):
        self.exs.append(ex)
        self.exs_str = str(ex) if self.exsLen() == 1 else self.exs_str + ',' + str(ex) 

    def addExs(self, exs):
        for ex in exs:
            if ex not in self.exs:
                self.addEx(ex)

    def getExs(self):
        return self.exs

    def exsLen(self):
        return len(self.exs)

    def exsKey(self):
        return self.exs_str

    def __str__(self):
        return "(%s):%d" % (self.exsKey(), self.pos)
    def __repr__(self):
        return self.__str__()

# A node that represents a segment in the segment graph
class SG_Node:
    key = SG_Key()
    width = 0
    end = 0

    isStart = False
    color = set()
    my_next = set()

    ntype = "E"

    def __init__(self, key, width, end):
        self.key = key
        self.width = width
        self.end = end
        self.my_next = set()
        self.color = set()
        self.ntype = "E"

    def addTxID(self, tx):
        self.color.add(tx)
       
    def getTxIDs(self):
        return(self.color)
        
    def keyStr(self):
        return self.key.__str__()
    
    def __repr__(self):
        return "%s --> width: %d, end:%d, ntype: %s, isStart: %s, color: (%s), next: %s" % \
              (self.keyStr(), self.width, self.end, self.ntype, self.isStart, ','.join(self.color), \
               self.my_next)

