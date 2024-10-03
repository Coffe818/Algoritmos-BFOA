
import random


def TumboNado2(secuencia):
    indices_guiones = [i for i, x in enumerate(secuencia) if x == '-']
    if not indices_guiones:
        return secuencia

    indice_guion = random.choice(indices_guiones)
    secuencia.pop(indice_guion)

    posiciones_validas = [i for i in range(len(secuencia) + 1)]

    nueva_posicion = random.choice(posiciones_validas)
    secuencia.insert(nueva_posicion, '-')

    return secuencia



def reproduccion(bacterias, puntuaciones):
    bacterias_puntuaciones = list(zip(bacterias, puntuaciones))
    bacterias_puntuaciones.sort(key=lambda x: x[1], reverse=True)
    bacterias_ordenadas, puntuaciones_ordenadas = zip(*bacterias_puntuaciones)

    n = len(bacterias)
    mitad = n // 2

    mejores_bacterias = list(bacterias_ordenadas[:mitad])
    mejores_puntuaciones = list(puntuaciones_ordenadas[:mitad])

    resultado_bacterias = mejores_bacterias 
    resultado_puntuaciones = mejores_puntuaciones



    return resultado_bacterias, resultado_puntuaciones


