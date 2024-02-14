import pyparsing as pp
from error_print import eprint
import global_vars
from functools import reduce

def multiline_parse_error(obj):
    eprint("Parsing error: Failed to parse block with label: {}".
           format(obj.label))
    exit(1)

data_pattern = pp.ZeroOrMore(pp.Group(pp.Word(pp.nums) +  pp.Suppress(':') + pp.common.number.add_parse_action(pp.common.convertToFloat)))

class Multiline:

    register = {}

    def __init__(self,label):
        self.label = label
        self.state = False
        self.dict = {}
        self.firstline = label + data_pattern
        __class__.register.setdefault(self.label, self)

        @self.firstline.set_parse_action
        def parseAction(s,loc,toks):
            self.lines.append(s.lstrip(self.label))
            self.enter()
            self.data = [None]*global_vars.ngtypes

        self.string = ''
        self.lines = []
        self.data_pattern = data_pattern.copy()
        
    def enter(self):
        assert self.state == False
        self.state = True

    def exit(self):
        assert self.state == True
        self.state = False

        @self.data_pattern.add_parse_action
        def dict_repr(toks):
            for k,v in toks:
                if k in self.dict:
                    multiline_parse_error(self)
                try:
                    self.dict[k] = v
                    self.data[int(k,global_vars.nalleles)] = v
                except ( KeyError, IndexError ):
                    multiline_parse_error(self)

        self.string += ' '.join(self.lines)
        self.data_pattern.parse_string(self.string)

        if None in self.data:
            multiline_parse_error(self)

multiline_labels = ['initial_state', 'fitness']

multiline_objects = [ Multiline(x) for x in multiline_labels ]

multiline_any = reduce(lambda x,y:x|y, [x.firstline for x in multiline_objects])

multiline_states = [ x.state for x in multiline_objects ]
