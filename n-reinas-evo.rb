#!/usr/bin/env ruby

# encoding: utf-8
# Programa: NReinasGenetico.rb
# Autor: Sebastián Narváez
# Email: sebasnr95@gmail.com
# Fecha creación: 2015-03-19
# Fecha última modificación: 2015-06-03
# Versión de Ruby: 2.2.2
# Versión del software: 0.3

# VERSIONES
# 0.1 - Se crea la clase NReinasGenetico y su constructor.
# 0.2 - Se crea la clase Cromosoma, con sus métodos mutar! y cruceUniforme
# 0.3 - Se crean los métodos evaluar población y seleccionarCromosomas de la 
#       clase NReinasGenetico

# Un cromosoma contiene un array de genes y uno de aptitud. Para el array de 
# genes, el index de una celda representa la fila del tablero, y su contenido 
# es un entero que representa # la columna en donde se encuentra la reina.
class Cromosoma

    attr_reader :cantReinas
    attr_accessor :genes, :aptitud, :id
    # 2, 0, 3, 1 ganador
    # Crea un cromosoma de tamaño cantReinas. Puede ser inicializado aleatoriamente
    # o crearse vacío. Retorna su id.
    def initialize(cantReinas, inicGenes = true)
        @genes = []
        @id = 0
        @aptitud = 0
        @cantReinas = cantReinas
        if(inicGenes)
            columnasOcupadas = []
            @cantReinas.times do |indice|
                columnaGenerada = Random.rand(@cantReinas)
                # Se verifica si la columna ya ha sido usada en este cromosoma,
                # para evitar cromosomas invalidos triviales.
                while columnasOcupadas[columnaGenerada] do
                    columnaGenerada = Random.rand(@cantReinas)
                end
                columnasOcupadas[columnaGenerada] = true
                @id += (2 ** indice) * columnaGenerada
                @genes.push columnaGenerada
            end
            #print @genes, " "
        end
    end

    # La mutación consiste en cambiar aleatoriamente dos posiciones del gen
    def self.mutar(cromosoma)
        nuevoCromosoma = Cromosoma.new(cromosoma.cantReinas, false)
        nuevoCromosoma.genes = cromosoma.genes.clone
        nuevoCromosoma.id = cromosoma.id
        pos1 = Random.rand(nuevoCromosoma.cantReinas)
        pos2 = Random.rand(nuevoCromosoma.cantReinas)
        if pos1 == pos2
            pos2 -= 1
        end
        # Se cambia el gen de la pos1 por el de la pos2
        nuevoCromosoma.genes[pos1], nuevoCromosoma.genes[pos2] = 
            nuevoCromosoma.genes[pos2], nuevoCromosoma.genes[pos1]
        # Se actualiza la id
        nuevoCromosoma.id -= (2 ** pos1) * cromosoma.genes[pos1]
        nuevoCromosoma.id -= (2 ** pos2) * cromosoma.genes[pos2]
        nuevoCromosoma.id += (2 ** pos1) * nuevoCromosoma.genes[pos1]
        nuevoCromosoma.id += (2 ** pos1) * nuevoCromosoma.genes[pos1]
        #print nuevoCromosoma.genes, " "
        return nuevoCromosoma
    end

    # El cruce uniforme genera individuos inválidos, por lo que no se usará
    def self.cruceUniforme(cromosomaPadre, cromosomaMadre)
        cromosomaHijo = Cromosoma.new(cromosomaPadre.cantReinas, false)
        porcentajeCruce = 0.6

        cromosomaHijo.cantReinas.times do |indexGen|
            indicador = Random.rand
            cromosomaHijo.genes.push ( indicador > porcentajeCruce ?
                                        cromosomaPadre.genes[indexGen] :
                                        cromosomaMadre.genes[indexGen] )
        end
		#print cromosomaHijo.genes, " - "
        return cromosomaHijo
    end
end

class NReinasGenetico

    # Se crean aleatoriamente @@tamPoblacion cromosomas.
    def initialize(cantReinas, cantGeneraciones = 50, tamPoblacion = 20, 
            evalCantAtaques = true, evalDiversidad = false, tamMatingPool = 2)
        #Inicializando el array de cromosomas y su tamaño
        @tamPoblacion = tamPoblacion
        @cantGeneraciones = cantGeneraciones
        @cantReinas = cantReinas
        @evalCantAtaques = evalCantAtaques
        @evalDiversidad = evalDiversidad
        @tamMatingPool = tamMatingPool

        # En un arreglo se guarda la población y en otro sus Ids. Los Ids son
        # calculados en el momendo de la creación de un cromosoma y ayudan 
        # en la verificación de la diversidad de la población.
        @poblacion = []
        @repetIds = Hash.new(0)
    end

    def ejecutar()
        matingPool = nil
        generacionFinal = 0;
        ataquesMejorCromosomaF = 0;
        @n_evaluaciones = 0
        @cantGeneraciones.times do |generacion|
            generacionFinal = generacion +1
            poblar(matingPool)
            @mejorCromosoma = @poblacion[0]
            #p "---------- GENERACION #{generacion} ----------"
            ataquesMejorCromosoma = self.calcularCantAtaques(@mejorCromosoma)
            ataquesMejorCromosomaF = ataquesMejorCromosoma
            if ataquesMejorCromosoma == 0
                #p "Se interrumpe el algoritmo en la generacion #{generacion}"
                break
            end
            matingPool = seleccionarCromosomas()
            #p "Repoblando..."
        end
        p "Generacion #{generacionFinal} - evaluacion #{@n_evaluaciones} - ataques #{ataquesMejorCromosomaF} - aptitud #{@mejorCromosoma.aptitud}"
        #p "#{@cantReinas}, #{@n_evaluaciones}, #{ataquesMejorCromosomaF}"
    end
    
    # Evalua la población y se calcula la aptitud de cada cromosoma respecto a:
    #   1 - Cantidad de ataques entre reinas (Por defecto).
    #   2 - Diversidad de los cromosomas.
    # Retorna el cromosoma con mejor aptitud de la población.
    def evaluarCromosoma(cromosoma)
		@n_evaluaciones += 1
        #mejorCromosoma = @poblacion[0]
        #p "Evaluando cromosoma: #{cromosoma.genes} respecto a:"
        cantAtaques = 0
        cantRepeticiones = 0

        if @evalCantAtaques
            #p "- Su cantidad de ataques"
            cantAtaques = self.calcularCantAtaques(cromosoma)
            #p "Hay #{cantAtaques} ataques"
        end
        if @evalDiversidad
            #p "- Su diversidad"
            # cantRepeticiones = calcularRepeticiones(cromosoma)
            cantRepeticiones = @repetIds[cromosoma.id]
            #p "Se encontraron #{cantRepeticiones} repticiones"
        end

        # Se calcula la aptitud del cromosoma
        cantAtaquesMax = ((@cantReinas - 1) * @cantReinas) / 2
        
        cromosoma.aptitud = 1 - ((Float(cantAtaques) / cantAtaquesMax) +
            (Float(cantRepeticiones) / @tamPoblacion))
            
        if cromosoma.aptitud > @mejorCromosoma.aptitud
            @mejorCromosoma = cromosoma
        end
        #p "Aptitud #{@mejorCromosoma.aptitud}"
        
        return @mejorCromosoma
    end

    def calcularCantAtaques(cromosoma)
		cantReinas = cromosoma.genes.length
        cantAtaques = 0
        # La reina de cada fila (excepto de la última)...
        (cantReinas - 1).times do |filaActual|
            # se contrasta con la reina de las demás filas y se verifica
            # si se cruzan
            (filaActual + 1).upto (cantReinas - 1) do |filaPorContrastar|
                # Para esto se haya la pendiente que forman las dos
                # reinas en el tablero. Si están en la misma diagonal,
                # la pendiente es 1 ó -1
                pendiente = Float(filaPorContrastar - filaActual) /
                            (cromosoma.genes[filaPorContrastar] -
                                cromosoma.genes[filaActual])
                if pendiente == 1 || pendiente == -1
                    #p "[#{filaActual}, #{cromosoma.genes[filaActual]}] y " +
                    #  "[#{filaPorContrastar}, #{cromosoma.genes[filaPorContrastar]}] se atacan"
                    cantAtaques += 1
                end
            end
        end

        return cantAtaques
    end

    # Realiza la seleccion de cromosomas de la poblacion actual por torneo. 
    def seleccionarCromosomas()
        #p "--- Empieza el torneo:"
        seleccion = []
        @tamMatingPool.times do
            competidor1 = @poblacion[Random.rand(@tamPoblacion)]
            competidor2 = @poblacion[Random.rand(@tamPoblacion)]
            evaluarCromosoma(competidor1)
            evaluarCromosoma(competidor2)
            ganador = competidor1.aptitud > competidor2.aptitud ? 
                        competidor1 : competidor2
            #p "--- Ganador: #{ganador.genes}"
            seleccion.push(ganador)
        end
        return seleccion
    end

    # Se vuelve a poblar a partir de un mating pool. Se realizan mutaciones
    # unicamente ya que los cruces producen 
    def poblar(matingPool = nil)
        @poblacion = []
        @tamPoblacion.times do
            if matingPool == nil
                nuevoCromosoma = Cromosoma.new(@cantReinas)
            else
                cromosomaAMutar = matingPool[Random.rand(@tamMatingPool)]
                #p "Empieza la mutacion de #{cromosomaAMutar.genes}"
                nuevoCromosoma = Cromosoma.mutar(cromosomaAMutar)
                #p "El resultado de la mutacion es #{cromosomaMutado.genes}"
            end
            @poblacion.push nuevoCromosoma
            @repetIds[nuevoCromosoma.id] += 1
        end
        #p ""
    end
end

## Aptitud por diversidad
#(3..30).each do |i|
#   NReinasGenetico.new(i, cantGeneraciones = 10000, tamPoblacion = 20, evalCantAtaques = false, evalDiversidad = true).ejecutar()
#end

#Aptitud por ataques
#(3..30).each do |i|
#   NReinasGenetico.new(i, cantGeneraciones = 1000, tamPoblacion = 20, evalCantAtaques = true, evalDiversidad = false).ejecutar()
#end

#Aptitud mixta
#(3..30).each do |i|
   NReinasGenetico.new(50, cantGeneraciones = 1000, tamPoblacion = 20, evalCantAtaques = true, evalDiversidad = true).ejecutar()
#end


