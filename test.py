t1 = ['Acidianus', 'Desulfurococcus', 'Ignicoccus']
t2 = ['Ignicoccus', 'Desulfurococcus', 'Acidianus']

missing_sp = list(set(t1)-set(t2))
print(missing_sp)
