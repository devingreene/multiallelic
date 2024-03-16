import pyparsing as pp
from error_print import eprint
import global_vars
from functools import reduce

TIMELIST_LENGTH = 50

def multiline_parse_error(obj):
    eprint("Parsing error: Failed to parse block with label: {}".
           format(obj.label))
    exit(1)

def k_correct(obj,k):
    assert isinstance(k,str) and k.isnumeric()
    if len(k) != global_vars.nloci or\
            any(not 0 <= int(c) < global_vars.nalleles for c in k):
                multiline_parse_error(obj)

data_patterns = {
        'dict': pp.ZeroOrMore(pp.Group(pp.Word(pp.nums) +  pp.Suppress(':') + pp.common.number.add_parse_action(pp.common.convertToFloat))),
        'list': pp.ZeroOrMore(pp.Word(pp.nums)),
        'timelist': pp.ZeroOrMore(pp.common.integer).add_parse_action(pp.common.convert_to_integer)
        }

class Multiline:

    register = {}

    def __init__(self,label,Type,default_value = None):
        self.label = label
        self.state = False
        self.default_value = default_value
        assert Type in ['dict', 'list', 'timelist']
        self.type = Type
        self.dict = {}
        self.set = set()
        self.firstline = label + data_patterns[self.type]
        __class__.register.setdefault(self.label, self)

        @self.firstline.set_parse_action
        def parseAction(s,loc,toks):
            # maxlen assigned here to give parser a chance to know what
            # ngtypes is.
            if self.type == 'list':
                self.maxlen = global_vars.ngtypes
            elif self.type == 'timelist':
                self.maxlen = TIMELIST_LENGTH
            self.lines.append(s.lstrip(self.label))
            multiline_labels_found.append(self.label)
            self.enter()
            if self.type == 'dict':
                self.data = [self.default_value]*global_vars.ngtypes
            elif self.type in ['list', 'timelist']:
                self.data = []
            else: assert False

        self.string = ''
        self.lines = []
        self.data_pattern = data_patterns[self.type].copy()
        
    def enter(self):
        assert self.state == False
        self.state = True

    def exit(self):
        assert self.state == True
        self.state = False

        if self.type == 'dict':
            @self.data_pattern.add_parse_action
            def dict_repr(toks):
                for k,v in toks:
                    k_correct(self,k)
                    if k in self.dict:
                        multiline_parse_error(self)
                    try:
                        self.dict[k] = v
                        self.data[int(k,global_vars.nalleles)] = v
                    except ( KeyError, IndexError ):
                        multiline_parse_error(self)
        elif self.type == 'list':
            @self.data_pattern.add_parse_action
            def list_repr(toks):
                for k in toks:
                    k_correct(self,k)
                    if k in self.set:
                        multiline_parse_error(self)
                    try:
                        self.set.add(k)
                        gtype = int(k,global_vars.nalleles)
                        if not 0 <= gtype < global_vars.ngtypes:
                            raise IndexError
                        self.data.append(gtype)
                    except ( KeyError, IndexError, ValueError ):
                        multiline_parse_error(self)
        elif self.type == 'timelist':
            @self.data_pattern.add_parse_action
            def timelist_repr(toks):
                if len(toks) > self.maxlen:
                    multiline_parse_error(self)
                self.data.extend(set(toks))
                self.data.sort()
        else: assert False

        self.string += ' '.join(self.lines)
        try:
            self.data_pattern.parse_string(self.string, parseAll=True)
        except pp.ParseException:
            multiline_parse_error(self)
        if self.type in ['list', 'timelist']:
            assert len(self.data) <= self.maxlen
            self.data.extend([-1] * (self.maxlen - len(self.data)))

        if None in self.data:
            multiline_parse_error(self)

# (label,Type,[default_value])
multiline_labels_types = [('initial_state','dict',0.),
                          ('fitness','dict',1.),
                          ('target_genotypes','list'),
                          ('report_points','timelist')]

multiline_objects = [ Multiline(*tup) for tup in multiline_labels_types ]

multiline_any = reduce(lambda x,y:x|y, [x.firstline for x in multiline_objects])

multiline_states = [ x.state for x in multiline_objects ]

multiline_labels_found = []
