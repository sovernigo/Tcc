import os

lista = os.listdir('testes_com_otimos/')
lista.sort()
#print(lista)
#exit()

for it in lista:
    subf = os.listdir(f'testes_com_otimos/{it}')
    for dados in subf:
        print(f'Executando {dados}')
        if '.out' not in dados:
            #print('oi')
            os.system(f'./exec testes_com_otimos/{it}/{dados} 100 100')