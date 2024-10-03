from bacteria import bacteria
from chemiotaxis import chemiotaxis
import numpy

poblacion = []
path = "../BFOA/multiFasta.fasta"
numeroDeBacterias = 5
numRandomBacteria = 1
iteraciones = 30
tumbo = 1                                              #numero de gaps a insertar 
nado = 3
chemio = chemiotaxis()
veryBest = bacteria(path)                #mejor bacteria   
tempBacteria = bacteria(path)            #bacteria temporal para validaciones
original = bacteria(path)                #bacteria original sin gaps
globalNFE = 0      #numero de evaluaciones de la funcion objetivo

dAttr= 0.0001 #0.1
wAttr= 0.0002 #0.2
hRep=dAttr
wRep= 10    #10



def clonaBest(veryBest, best):
    veryBest.matrix.seqs = numpy.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction
    

def validaSecuencias(path, veryBest):
    #clona a veryBest en tempBacteria   
    tempBacteria.matrix.seqs = numpy.array(veryBest.matrix.seqs)
    #descartar los gaps de cada secuencia
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-","")
    #tempBacteria.tumboNado(1)    

    #valida que las secuencias originales sean iguales a las secuencias de tempBacteria
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return

def tumboNado_iterativo(bacteria, tumbo, pasos=5):
    #print(f"Iniciando quimiotaxis para la bacteria con fitness inicial: {bacteria.fitness}")
    
    for paso in range(pasos):
        #print(f"  Paso  {paso + 1}/{pasos}")
        
        # 1. Aplicar tumbo/nado a la bacteria
        nueva_bacteria = bacteria.clonar(path)  # Clonar la bacteria original para modificarla
        nueva_bacteria.tumboNado(tumbo)
        nueva_bacteria.autoEvalua()  # Evaluar la nueva configuración
        
        # 2. Comparar la nueva puntuación con la original
        if nueva_bacteria.fitness > bacteria.fitness:
            # Si mejora, actualizar la bacteria original
            #print(f"  Mejora encontrada con fitness: {nueva_bacteria.fitness}")
            bacteria.matrix.seqs = nueva_bacteria.matrix.seqs
            bacteria.fitness = nueva_bacteria.fitness
        # else:
        #     print(f"  No hay mejora, manteniendo fitness: {bacteria.fitness}")
    return bacteria  # Retornar la bacteria mejorada


for i in range(numeroDeBacterias):                                            #poblacion inicial
    poblacion.append(bacteria(path))
    #print("bacteria: ",i,poblacion[i].matrix.seqs)


for _ in range(iteraciones  ):                                                  #numero de iteraciones  
    for i, bacteria_actual in enumerate(poblacion):
        poblacion[i] = tumboNado_iterativo(bacteria_actual, tumbo)

        #print("blosumScore: ",bacteria.blosumScore)
    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)                 #d_attr, w_attr, h_rep, w_rep):
    globalNFE += chemio.parcialNFE 
    best = max(poblacion, key=lambda x: x.fitness)
    if (veryBest == None) or (best.fitness > veryBest.fitness):
        clonaBest(veryBest, best)
    print("interaccion: ",veryBest.interaction,"fitness: ",veryBest.fitness, " NFE:",globalNFE )
    #veryBest.showGenome()
    chemio.eliminarClonar(path, poblacion)
    chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)                #inserta  bacterias aleatorias
    print("poblacion: ",len(poblacion))

print("Mejor bacteria: ")
veryBest.showGenome()
validaSecuencias(path, veryBest)