import Ref.TP53Ref as TP53Ref
import Ref.GATA3Ref as GATA3Ref
import Ref.PIK3CARef as PIK3CARef
import Ref.CDH1Ref as CDH1Ref
import random

var = "T"
start = 68857481

s = CDH1Ref.CDH1Ref
i = start - 68771195
print(s[i])
s = s[:i] + var + s[i + 1 :]
print(s[i - random.randint(1, 100) : i + random.randint(1, 100)])
