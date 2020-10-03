# Metropolis-Ising2D
Algoritmo para calcular energía y magnetización de una red de espines segun el modelo de Isign

Este algoritmo se basa en el método de Metropolis de aceptación o rechazo de estados generados aleatoreamente pesados por una distribución de probabilidad tipo Boltzman. El agoritmo fue escrito en lenguaje C y compilado bajo Linux-Ubuntu 18.04 con la orden g++ -o Ising.exe Ising2D.C Functions.C Scripts.C Files.C; time ./Ising.exe la cual corre el código y llama las funciones necesarias para su uso. Adicionalmente se necesita tener instalado GNUPlot para el procesado de gráficos pues el código llama al programa para generar los gráficos de energia y magnetización. El archivo scritp permite modificar los parámetros de GNUPlot para el tipo de gráficos.
