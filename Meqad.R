#======================================================================================
# Author: Franze Costa
# Universidade Federal da Paraiba
#======================================================================================
#     INFORMAcoES PRELIMINARES
#======================================================================================
Meqad=function(){
	print("Esse arquivo contem as rotinas para diversos resultados e simulacoes.")
	print("Em geral, a maior parte das funcoes ja esta implementada no R.")
	print("O que temos aqui e uma organizacao alternativa.")
	print("Para  maiores detalhes, basta solicitar detalhes pela funcao 'mqd()")
	print("As funcoes sao organizadas pelos numeros abaixo indicados:")
	print("1 - Para geracao de graficos e medidas descritivas")
	print("2 - Para analise de quantis e distribuicoes")
	print("3 - Para analise de correlacao")
	print("4 - Para manuseio de outliers e missing values")
	print("5 - Rotinas para modelagem linear")
	print("6 - Rotinas para analise de variância")
	print("7 - Simulacao de estimadores de algumas distribuicoes de probabilidades")
}

mqd=function(fun=1){ 
	if(fun==1){
		print("Temos abaixo instrucoes de uso das rotinas Meqad de analise descritiva")
		print("-------------------------------------------------------------------------------")
		print("Os comandos sao:")

		print("apara(DADOS)         Apara os dados dos valores extremos (default de 10%)")
		print("winsoriza(DADOS)     Winsoriza os dados de uma variavel (default de 10%)")
		print("graf.mqd(DADOS)      Gera quatro graficos de verificacao exploratoria")
		print("posi.mqd(DADOS)      Apresenta medidas de posicao")
		print("disp.mqd(DADOS)      Apresenta medidas de dispersao")
		print("assi.mqd(DADOS)      Apresenta medidas de assimetria")
		print("curt.mqd(DADOS)      Apresenta medidas de curtose")
		print("summary.mqd(DADOS)   Extrai todo o conjunto de graficos e medidas")
	}
	else if(fun==2) {
		print("Temos abaixo instrucoes de uso das rotinas Meqad de analise de quantis")
		print("-------------------------------------------------------------------------------")
		print("Os comandos sao:")
		
		print("COMANDOS PARA INTERVALOS DE CONFIANcA DE QUANTIS")
		print("ic.q.mqd(x, q)           Indica a estimativa e o intervalo de confianca do quantil q")
		print("ic.decis.mqd(x)          Estimativas e o intervalos de confianca dos nove decis")
		print("ic.central.mqd(x)        Estimativas e o intervalos de confianca para posicao central")
		
		print("COMANDOS PARA COMPARAcaO DE QUANTIS DE DUAS VARIaVEIS OU AMOSTRAS")
		print("gq.mqd(x,y)              Gera medidas descritivas e o grafico de decis ")
		print("quniv.test.mqd(x, med, q)Teste (exato) para o valor de um quantil q de UMA variavel")
		print("central.test.mqd(x)      Testes t, de Wilcoxon e da mediana para a posicao central")
		print("mediana.test.mqd(x,y)    Testes para comparacao de medianas de duas amostras")
		print("qbiv.test.mqd(x, y, q)   Gera o teste para um quantil especifico q de duas amostras")	
		print("decis.test.mqd(x, y)     Gera o teste para os nove decis e o grafico")
		print("q.test.mqd(x, y, med, q) Gera geral para uma ou duas amostras - quinv e qbiv")
		print("cvm.mqd(x, y)            Realiza o teste de Cramer-von Mises")
		print("dist.test.mqd(x, y)      Realiza os testes de igualdade de duas distribuicoes")
		print("gera.quantis.mqd(x,y)    Gera quantis de k variaveis e o grafico conjunto")
		 
		print("--------------------------------------------------------------------------------------")
		print("COMPLEMENTO - ANOVA DE QUANTIS")
		
		print("Suposicao - x e y sao pareadas")
		print("1 - Tomar uma variavel quantitativa x")
		print("2 - Tomar uma variavel categorica y")
		print("3 - Conferir se x e quantitativa - is.numeric(x)")
		print("4 - Conferir se y e categorica - is.factor(y)")
		print("5 - Definir o o data-frame - z=data.frame(x,y)")
		print("6 - Definir a variavel de comparacao - z=fac2list(z$x,z$y)")
		print("7 - Implementar a anova quantilica para um quantil q - Qanova(z, q)")
		print("8 - Proceder a descricao: ")
		print("  8.1. Isolar as duas variaveis z1=x[y==...]; z2=x[y==...]; ...")
		print("  8.2. Comparar os intervalos de confiana: ic.q.mqd(z1, q); ic.q.mqd(z2, q); ...")
	}
	else if(fun==3) {
		print("Temos abaixo instrucoes de uso das rotinas Meqad de Analise de correlacao")
		print("-------------------------------------------------------------------------------")
		print("Os comandos sao:")
		print("cor.w.mqd            Gera correlacao winsorizada")
		print("scor.mqd             Gera correlacao com exclusao de outliers (skipped correlation)")
		print("cor.mqd(x, y)        Gera medidas de correlacao e o p-valor")
	}
	else if(fun==4) {
		print("Temos abaixo instrucoes de uso das rotinas Meqad de simulacao de ICs")
		print("-------------------------------------------------------------------------------")
		print("Os comandos sao:")
		print("exc.mv.mqd(x)          Exclui as linhas de um base de dacis que possuen missing values")
		print("imp.mv.mqd(x)          Imputa missing values (media ou estocastico) em uma variavel especifica")
		print("exc.out.mqd(x)         Exclui outliers. Criterios: interquartil (Carling), convencional e mad-madiana")
	}
	else if(fun==5) {
		print("Temos abaixo instrucoes de uso das funcoes de apoio para modelagem linear")
		print("-------------------------------------------------------------------------------")
		print("Os comandos sao:")
		print("jarque.bera.mqd(x)           Para o teste de normalidade de Jarque-Bera")
		print("normalidade.mqd(x)           Para os principais teste de normalidade")
		print("homoscedasticidade.mqd(x)    Para os principais teste de homoscedasticidade")
		print("independencia.mqd(x)         Para os principais teste de independência")
		print("diagplot.mqd(modelo)         Para quatro graficos de diagnostico")
		print("env.mqd(x)                   Para o envelope simulado")
		print("mnl.mqd(formula)             Modelagem normal linear completa")
		print("confint.rfit.mqd(modelo)     Intervalo de confianca do modelo rfit")
		print("mbr.mqd(formula)             Modelagem por ranques completa")
		print("mrq.mqd(formula)             Modelagem quantilica completa")
	}
	else if(fun==6) {
		print("geral.oneway.mqd(resp, categ)      Gera os resultados completos dos três testes")
		print("posthoc.oneway.mqd(resp, categ)    Gera os resultados da comparacao dois a dois")
		print("descritivo.oneway.mqd(resp, categ) Gera as medidas descritivas pelas categorias")
		print("oneway.mqd(resp, categ)            Desenvolve anova parametrica e nao parametrica")
	}
	else if(fun==7) {
		print("Temos abaixo instrucoes de uso das simulacoes implementadas")
		print("-------------------------------------------------------------------------------")
		print("Os comandos sao:")
		print("simula.bern.mqd(p)           Simula intervalos de confianca para o parâmetro da distribuicao de Bernoulli")
		print("simula.unif.mqd(min, max)    Simula intervalos de confianca para os parâmetros da distribuicao uniforme")
		print("simula.pois.mqd(lambda)      Simula intervalos de confianca para o parâmetro da distribuicao Poisson")
		print("simula.rexp.mqd(lambda)      Simula intervalos de confianca para o parâmetro da distribuicao exponencial")
		print("simula.norm.mqd(mu, sigma2)  Simula intervalos de confianca para os parâmetros da distribuicao normal")
	}
}

options( warn = -1 )
#======================================================================================
#     ROTINA DE ANaLISE DESCRITIVA
#======================================================================================
library(ggplot2)
#--------------------------------------------------------------------------------------
#EXPLORAcaO INICIAL
#--------------------------------------------------------------------------------------
graf.mqd=function(x) {
	print(paste("GRaFICOS: boxplot, histograma e quantis"))
	par(mfrow=c(2,2)) 
  	boxplot(x, ylab="Valores"); title("Boxplot")							
  	hist(x, main="Histograma da variavel", xlab="Valores", ylab="Frequência");							
  	qqnorm(x, main="QQPlot da variavel", xlab="Quantis da normal", ylab="Quantis observados")
	qqline(x)		
  	plot(seq(0:10), quantile(x, (0:10)/10), type="b", xlab="Decis", ylab="Valores observados")
	title("Grafico de decis")
}

#--------------------------------------------------------------------------------------
#ALGORITMOS PRELIMINARES
#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
#MEDIDA DE BIWEIGHT MIDVARIANCE
mw.var.mqd=function(x){
	n=length(x)
	M=median(x)
	u=(x-M)/(9*mad(x))
	ai=ifelse(abs(u)<1,1,0)	
	num=n*sum(ai*(x-M)^2*(1-u^2)^4)
	den=sum(ai*(1-u^2)*(1-5*u^2))^2
	var=num/den
	var
}

#--------------------------------------------------------------------------------------
#PROCESSO DE APARA
apara=function(x, a=0.1) {
	aux=sort(x); n=length(aux)
	if (a<=0.5) {
		pinf=floor(a*n)+1; qinf=aux[pinf]
		psup=n-pinf+1; qsup=aux[psup]
		x=ifelse(x<qinf, NA, x)
		x=ifelse(x>qsup, NA, x)
		x=data.frame(x)
	}
	else print("ERRO!! O nivel maior que 0.5 nao permite aparas") 
	x=exc.mv.mqd(x)
	x
}

#--------------------------------------------------------------------------------------
#PROCESSO DE WINSORIZAcaO 
winsoriza=function(x, w=0.1) {
	aux=sort(x); n=length(aux)
	if (w<0.5) {
		pinf=floor(w*n)+1; qinf=aux[pinf]
		psup=n-pinf+1; qsup=aux[psup]
		x=ifelse(x<qinf, qinf, x)
		x=ifelse(x>qsup, qsup, x)
		x
	}
	else print("ERRO!! O nivel maior que 0.5 nao permite winsorizacao") 
}

#-----------------------------------------------------------------------
#MeDIA WINSORIZADA
media.w.mqd=function(x, w=0.1){
	y=winsoriza(x, w)
	mean(y)
}

#-----------------------------------------------------------------------
#VARIÂNCIA E DESVIO APARADOS E WINSORIZADOS
#Desvio aparado
var.a.mqd=function(x, a=0.1){
	y<-apara(x)
	var(y)
}

desv.a.mqd=function(x, a=0.1){
	y<-apara(x)
	sd(y)
}

var.w.mqd=function(x, w=0.1){
	y=winsoriza(x, w)
	var(y)
}

desv.w.mqd=function(x, w=0.1){
	y=winsoriza(x, w)
	sd(y)
}

#--------------------------------------------------------------------------------------
#MEDIDAS DE POSIcaO
#--------------------------------------------------------------------------------------
posi.mqd=function(x) {
	print(paste("MEDIDAS DE POSIcaO"))
	m_simples=mean(x)
	n=length(x)
	des=sd(x)
	lim=m_simples+qt(.025, n-1)*des/sqrt(n)
	lsm=m_simples+qt(.975, n-1)*des/sqrt(n)
	m_apa_10=mean(x, 0.1)
	Min=min(x)
	Max=max(x)
	quart_1=quantile(x, .25)
	quart_2=quantile(x, .50)
	quart_3=quantile(x, .75)
	win_1=media.w.mqd(x, 0.1)

	print(paste("1 - Media aritmetica simples:", round(m_simples, 3)))
	print(paste("2 - Erro padrao da media:", round(des/sqrt(n), 3)))
	print(paste("3 - Intervalo de confianca(padrao) de 95%:", "LI-", round(lim, 3), "; LS-", round(lsm, 3)))
	print(paste("4 - Media aparada de 10%:", round(m_apa_10, 3)))
	print(paste("5 - Media winsorizada de 10%:", round(win_1, 3)))
	print(paste("6 - Minimo e maximo:", "Min-", round(Min, 3), "; Max-", round(Max, 3)))
	print(paste("7 - Quartis:", "Q1-", round(quart_1, 3),"; Q2-", round(quart_2, 3), "; Q3-",round(quart_3, 3)))
}

#--------------------------------------------------------------------------------------
#MEDIDAS DE DISPERSaO
#--------------------------------------------------------------------------------------
disp.mqd=function(x) {
	print(paste("MEDIDAS DE DISPERSaO"))
	at=max(x)-min(x)
	aiq=quantile(x,0.75)-quantile(x,0.25)
	saiq=aiq/2
	dm=mean(abs(x-mean(x)))
	v=var(x)
	dp=sd(x)
	disp_trim=function(x, q) {
		novo_x=x[x>=quantile(x, q) & x<=quantile(x,1-q)]
		sd(novo_x)
	}
	va10=disp_trim(x, .1)^2
	vw1=var.w.mqd(x, 0.1)
	bi.mid=mw.var.mqd(x)
	dma=mad(x)
	print(paste("1 - Amplitude total:", round(at, 3)))
	print(paste("2 - Intervalo interquartil semi interquartil:", "IQ-", round(aiq, 3), "; SIQ-",round(saiq, 3)))
	print(paste("3 - Desvio medio:", round(dm, 3)))
	print(paste("4 - Variância e desvio padrao:", round(v, 3), ";", round(v^.5, 3)))
	print(paste("5 - Var. e desvio aparados de 10%:", round(va10, 3), ";", round(va10^.5, 3)))
	print(paste("6 - Var. e desvio winsorizados de 10%:", round(vw1, 3), ";", round(vw1^.5, 3)))
	print(paste("7 - Biweight midvariance e desvio:", round(bi.mid, 3), ";", round(bi.mid^.5, 3)))
	print(paste("8 - Desvio mediano absoluto - MAD:", round(dma, 3)))
}

#--------------------------------------------------------------------------------------
#MEDIDAS DE ASSIMETRIA				
#--------------------------------------------------------------------------------------
assi.mqd=function(x) {
	print(paste("COEFICIENTES DE ASSIMETRIA"))
	sk1=(sum((x-mean(x))^3)/length(x))/sd(x)^3
	sk2=((quantile(x,0.75)-median(x))-(median(x)-quantile(x,0.25)))/(quantile(x,0.75)-quantile(x, 0.25))
	sk3=((quantile(x,0.90)-median(x))-(median(x)-quantile(x,0.10)))/(quantile(x,0.90)-quantile(x, 0.10))
	sk_80=((quantile(x,0.80)-median(x))-(median(x)-quantile(x,0.20)))/(quantile(x,0.80)-quantile(x, 0.20))
	sk_95=((quantile(x,0.95)+quantile(x,1-0.95)-2*median(x)))/(quantile(x,0.95)-quantile(x,1-0.95))

	print(paste("Obs: em todas as medidas, zero indica simetria"))
	print(paste("1 - Assimetria de Pearson:", round(sk1, 3)))
	print(paste("2 - Assimetria de Bowley:", round(sk2, 3)))
	print(paste("3 - Assimetria de Kelley:", round(sk3, 3)))
	print(paste("4 - Coef. quantilico de 0.80:", round(sk_80, 3)))
	print(paste("5 - Coef. quantilico de 0.95:", round(sk_95, 3)))
}

#--------------------------------------------------------------------------------------
#MEDIDAS DE CURTOSE			
#--------------------------------------------------------------------------------------
#COEFICIENTE MOMENTO DE CURTOSE
curt.mqd=function(x) {
	print(paste("COEFICIENTES DE CURTOSE"))
	kurt1=(sum((x-mean(x))^4)/length(x))/sd(x)^4
	kurt2=(quantile(x, 0.75)-quantile(x, 0.25))/(2*(quantile(x, 0.90)-quantile(x,0.10)))

	print(paste("1 - Curtose de Pearson (3 e o valor da distribuicao mesocurtica):", round(kurt1, 3)))
	print(paste("2 - Coefic. quantilico (0,263 e o valor da distribuicao mesocurtica):", round(kurt2, 3)))
}

#--------------------------------------------------------------------------------------
#VERIFICAcaO CONJUNTA			
#--------------------------------------------------------------------------------------
#COEFICIENTE MOMENTO DE CURTOSE
summary.mqd=function(DADOS) {
	graf.mqd(DADOS)		#Gera quatro graficos de verificacao exploratoria
	print(paste("---------------------------------------------------"))
	posi.mqd(DADOS)		#Apresenta medidas de posicao
	print(paste("---------------------------------------------------"))
	disp.mqd(DADOS)		#Apresenta medidas de dispersao
	print(paste("---------------------------------------------------"))
	assi.mqd(DADOS)		#Apresenta medidas de assimetria
	print(paste("---------------------------------------------------"))
	curt.mqd(DADOS)		#Apresenta medidas de curtose
}


#======================================================================================
#    ANaLISE DE OUTLIERS E MISSING VALUES
#======================================================================================
#-----------------------------------------------------------------------
#REMOcaO DE OUTLIERS
#-----------------------------------------------------------------------
exc.out.mqd=function(x, method="interquartil") {
#Exclui outliers de uma variavel x, transformando em dado perdido
#Metodo interquartil: outliers e o valor que estiver fora do intervalo definido pela mediana
# mais ou menos k vezes o intervalo interquartil. k e a o valor definido por Carling e e funcao 
# do tamanho da amostra. Referência: Carling, K. Resistant outlier rules and the non-Gaussian case. 
# Computational Statistics & Data Analysis, 33(3), 249-258, 2000.

#Metodo convencional: outliers e o valor que estiver fora do intervalo definido pela media mais ou menos 
# 2.24 vezes o desvio padrao

#Metodo mad-mediana: outliers e o valor que estiver fora do intervalo definido pela mediana mais ou menos 
# 2.24 vezes o desvio mediano absoluto

	if(method=="interquartil") {
		n=length(x)
		k=(17.63*n-23.64)/(7.74*n-3.71)
		li=median(x)-k*(quantile(x, .75)-quantile(x, .25))
		ls=median(x)+k*(quantile(x, .75)-quantile(x, .25))
		rem1<-c(x[x<li])
		rem2<-c(x[x>ls])
		remove=c(rem1, rem2)
		x[is.element(x, remove)]<-NA
	}
	else if(method=="convencional") {
		li=mean(x)-2.24*sd(x)
		ls=mean(x)+2.24*sd(x)
		rem1<-c(x[x<li])
		rem2<-c(x[x>ls])
		remove=c(rem1, rem2)
		x[is.element(x, remove)]<-NA
	
	}
	else if(method=="mad-mediana") {
		li=median(x)-2.24*(mad(x)/.6745)
		ls=median(x)+2.24*(mad(x)/.6745)
		rem1<-c(x[x<li])
		rem2<-c(x[x>ls])
		remove=c(rem1, rem2)
		x[is.element(x, remove)]<-NA
	}
	x
}

#-----------------------------------------------------------------------
#REMOcaO DE MISSING VALUES
#-----------------------------------------------------------------------
exc.mv.mqd=function(dados) {
#Remove valores perdidos das linhas de uma base de dados em que ha ao menos 
# um dado perdido

	n1=nrow(dados)
	coluna_aux=c(1:nrow(dados))
	for(i in 1:nrow(dados)) {
		if(sum(is.na(dados[i,])>=1)) coluna_aux[i]<-0
	}
	sem_missing=dados[coluna_aux[coluna_aux>0],] 
	sem_missing
}

#-----------------------------------------------------------------------
#IMPUTAcaO DE VALORES A MISSING VALUES
#-----------------------------------------------------------------------
imp.mv.mqd=function(x, method="media") {
#Imputacao de valores a missing values de uma dada variavel
#Metodos: media da variavel e estocastico (a media mais um valor aleatorio
#  extraido de uma distribuicao normal com media zero e a variância amostral

	aux=c(1:length(x))
	media=sd(x, na.rm=TRUE)
	desvio=sd(x, na.rm=TRUE)
		if(method=="media") {
			for (i in 1:length(x)) {
				if(is.na(x[i])==TRUE) aux[i]<-media
				else aux[i]<-x[i]
			}
		}
		if(method=="estocastico") {
			for (i in 1:length(x)) {
				if(is.na(x[i])==TRUE) aux[i]<-media+rnorm(1,0,desvio)
				else aux[i]<-x[i]
			}	
		}
	x=aux
	x
}

#======================================================================================
#    ANaLISE DE CORRELAcaO
#======================================================================================

#-----------------------------------------------------------------------
#CORRELAcaO WINSORIZADA
#-----------------------------------------------------------------------
cor.w.mqd=function(x, y, w=0.1, sig=0.05) {
	n=length(x)
	mat=cbind(x,y)
	x1=winsoriza(mat[,1], w)
	y1=winsoriza(mat[,2], w)
	mat1=cbind(x1,y1)
	rw=cor(x1, y1)
	tw=rw*((n-2)/(1-rw^2))^0.5
	glw=n-floor(n*w)-2
	sigw=2*(min(pt(tw, glw), 1-pt(tw, glw)))
	TH=ifelse(sigw<sig, "rejeitamos a hipotese de correlacao nula",
		"nao rejeitamos a hipotese de correlacao nula")
	print(paste("A correlacao winsorizada de", 100*w, "% e", round(rw, 3)))
	print(paste("O p-valor e:", round(sigw, 4)))
	print(paste("A", 100*sig, "%,", TH))
}


#-----------------------------------------------------------------------
#CORRELAcaO COM REMOcaO DE OUTLIER (SKIPPED)
#-----------------------------------------------------------------------
scor.mqd=function(x,y, method=c("interquartil", "convencional", "mad-mediana")) {
#Metodo default: mediana e intervalo interquartil.
	n=length(x)
	x_sem_out<-exc.out.mqd(x)
	y_sem_out<-exc.out.mqd(y)
	scor=cor(x_sem_out, y_sem_out, "complete.obs")
	par(mfrow=c(3,2))
	hist(x, na.rm=TRUE); hist(x_sem_out, na.rm=TRUE)
	hist(y, na.rm=TRUE); hist(y_sem_out, na.rm=TRUE)
	plot(x,y, na.rm=TRUE); plot(x_sem_out, y_sem_out, na.rm=TRUE)
	res=list(Pearson=cor(x,y), Skipped=scor, elim_x=n-length(which(!is.na(x_sem_out))), elim_y=n-length(which(!is.na(y_sem_out))))
	res
}


#-----------------------------------------------------------------------
#DIVERSAS CORRELAcoES
#-----------------------------------------------------------------------
cor.mqd=function(x,y, w=.1) {
	print(paste("GERA MEDIDAS DE CORRELAcaO: WINSORIZADA, PEARSON, SPEARMAN E KENDALL"))	
	n=length(x)	
	#PROCEDER a WINSORIZAcaO
	winsoriza=function(x, ww=w) {
		aux=sort(x)
		if (w<0.5) {
			pinf=floor(w*n)+1; qinf=aux[pinf]
			psup=n-pinf+1; qsup=aux[psup]
			x=ifelse(x<qinf, qinf, x)
			x=ifelse(x>qsup, qsup, x)
			x
   		}
		else print("ERRO!! O nivel maior que 0.5 nao permite winsorizacao") 
	}
	#CALCULANDO A CORRELAcaO WINSORIZADA
	cor.w.mqd=function(x, y, wc=w, sig=0.05) {
		mat=cbind(x,y)
		x1=winsoriza(mat[,1], wc)
		y1=winsoriza(mat[,2], wc)
		mat1=cbind(x1,y1)
		rw=cor(x1, y1)
		tw=rw*((n-2)/(1-rw^2))^0.5
		glw=n-floor(n*w)-2
		sigw=2*(min(pt(tw, glw), 1-pt(tw, glw)))
		med=round(cbind(rw, sigw),3)
		print(paste("Cor. winsor. de", 100*w, "% (default - 10%) - valor:", med[,1], "; p-valor:", med[,2]))
	}	

	#CALCULANDO A CORRELAcaO SEM OUTLIERS
	x_sem_out<-exc.out.mqd(x)
	y_sem_out<-exc.out.mqd(y)
	scor=cor.test(x_sem_out, y_sem_out, na.rm=TRUE)
	skcor=round(scor$estimate, 3)
	skpval=round(scor$p.value, 3)
	scor$parameter
	print(paste("Skipped correlation - valor:", skcor, "; p-valor:", skpval, "; pares eliminados:", n-scor$parameter))
	
	#EXTRAINDO CORRELAcoES CONVENCIONAIS
	cp=cor.test(x,y, method="pearson"); 
	corp=round(cp$estimate, 3)
	pvalp=round(cp$p.value, 3)
	cs=cor.test(x,y, method="spearman")
	cors=round(cs$estimate, 3)
	pvals=round(cs$p.value, 3)
	ck=cor.test(x,y, method="kendall")
	cs=cor.test(x,y, method="spearman")
	cork=round(ck$estimate, 3)
	pvalk=round(ck$p.value, 3)
	cw=cor.w.mqd(x,y)
	print(paste("Cor. de Pearson - valor:", corp, "; p-valor:", pvalp))
	print(paste("Cor. de Spearman - valor", cors, "; p-valor:", pvals))
	print(paste("Cor. de Kendall - valor", cork, "; p-valor:", pvalk))
}

#======================================================================================
#     ROTINA DE INTERVALO DE CONFIANcA DE QUANTIS
#======================================================================================

#--------------------------------------------------------------------------------------
#QUANTIL DE HARRELL-DAVIS
#--------------------------------------------------------------------------------------
hdquantil.mqd=function(x, q=.5) {
	if(q>1) stop("O quantil tem que ser definido entre 0 e 1")
	if(q<0) stop("O quantil tem que ser definido entre 0 e 1")
	x<-sort(x)
	n=length(x)	
	n1=n+1
	hdquantil=c()
	for(i in 1:length(q)){
		wi=pbeta(1:n/n, q[i]*n1, (1-q[i])*n1)-pbeta(0:(n-1)/n, q[i]*n1, (1-q[i])*n1)
		hdquantil[i]=sum(wi*x)
	}
	round(hdquantil[], 4)
}

#--------------------------------------------------------------------------------------
#IC PARA UM QUANTIL QUALQUER
#--------------------------------------------------------------------------------------
ic.q.mqd=function(x, q=0.5, b=0, tamanho=95) {
#Para n ate 20 e usado o intervalo de confianca exato
#Para n maior que 20 e o usado o intervalo aproximado da binomial pela normal
#Opcao de boostrapping se b>0. Sugestao de acima de 500.
	xo=sort(x)
	n=length(x)
	alfa=(100-tamanho)/100
	z=qnorm(alfa/2)
	Quantil<-hdquantil.mqd(x, q)
	if(b>0) {
		quant=c()
	  	for(i in 1:b) {
	  		quant[i]=hdquantil.mqd(sample(x, n, replace=TRUE), q)
	  	}
	  	ic=round(c(Quantil=Quantil, LI=hdquantil.mqd(quant, .025), LS=hdquantil.mqd(quant, .975)), 3)
	}
	else {
		if(n<=20) {
			p.inf=max(1, qbinom(.025, n, q))
			p.sup=min(n, qbinom(.975, n, q))
			ic=round(c(Quantil=Quantil, LI=xo[p.inf], LS=xo[p.sup]), 3)
		}
		else {
			p.inf=max(1, floor(n*q+.5+z*sqrt(n*q*(1-q))))
			p.sup=min(n, ceiling(n*q+.5-z*sqrt(n*q*(1-q))))
			ic=round(c(Quantil=Quantil, LI=xo[p.inf], LS=xo[p.sup]), 3)
		}
	}
	ic
}

#--------------------------------------------------------------------------------------
#PARA OS DECIS
#--------------------------------------------------------------------------------------
ic.decis.mqd=function(x, b=0) {
	dec=matrix(1,9, 4)
	for (i in 1:9) {
		dec[i,1]=i/10
		dec[i, 2:4]=ic.q.mqd(x, i/10, b=b)
	}

	Resq=as.data.frame(round(dec, 3))
	names(Resq)=c("Decil", "Valor", "LI", "LS") 
	plot(Resq[,2], main="Grafico de decis", xlab="Decis", ylab="Valores", ylim=c(Resq[1,3],Resq[9,4]),
		type="c"); text(Resq[,2], label=(Resq[,2]), col="red", cex=0.7)
	lines(Resq[,3], type="c", col="red"); text(Resq[,3], label=(Resq[,3]), cex=0.7)
	lines(Resq[,4], type="c", col="red"); text(Resq[,4], label=(Resq[,4]), cex=0.7)
	Resq
}

#--------------------------------------------------------------------------------------
#PARA OS MeDIA E MEDIANA - POR BOOSTRAPPING
#--------------------------------------------------------------------------------------
ic.central.mqd<-function(x, b=1000, tamanho=95) {
	print(paste("Extracao por boostrap com", b, "iteracoes, e nivel de confianca de", tamanho, "%"))
	print("Nivel de apara e winsorizacao de 10%")
	
	alfa=(100-tamanho)/100
	n=length(x)
	media.apa=c()
	media.win=c()
	media=c()
	mediana=c()
	for (i in 1:b) {
		media[i]<-mean(sample(x, n, replace=TRUE))
		media.apa[i]<-mean(sample(x, n, replace=TRUE), 0.1)
		media.win[i]<-media.w.mqd(sample(x, n, replace=TRUE), 0.1)
		mediana[i]<-quantile(sample(x, n, replace=TRUE), .5)
	}
	IC_media=c(Valor=mean(media), LI=quantile(media, alfa/2), LS=quantile(media, 1-alfa/2))
	IC_media.apa=c(Valor=mean(media.apa), LI=quantile(media.apa, alfa/2), LS=quantile(media.apa, 1-alfa/2))
	IC_media.win=c(Valor=mean(media.win), LI=quantile(media.win, alfa/2), LS=quantile(media.win, 1-alfa/2))
	IC_mediana=c(Valor=mean(mediana), LI=quantile(mediana, alfa/2), LS=quantile(mediana, 1-alfa/2))
	resultado=data.frame(cbind(round(IC_media, 3), round(IC_media.apa, 3), round(IC_media.win, 3), round(IC_mediana, 3)))
	names(resultado)<-c("Med. simples", "Med. aparada", "Med. winsor.", "Mediana")
	resultado
}

#======================================================================================
#     ROTINA DE COMPARAcaO DE QUANTIS
#======================================================================================
#---------------------------------------------------------------------------
#PARTE 1 - ANaLISE DESCRITIVA GRaFICA
#---------------------------------------------------------------------------
gq.mqd=function(x,y) {
#Gera grafico de decis de DUAS amostras, independentes ou nao
	q=(1:9)/10	
	q_x=hdquantil.mqd(x, q)
	q_y=hdquantil.mqd(y, q)
	mi=min(c(x, y))	
	ma=max(c(x, y))
	val=as.data.frame(round(cbind(q_x, q_y), 3))
	names(val)=c("Decis da amostra 1", "Decis da amostra 2") 

	#Visualizacao grafica
	plot(q, q_x, main="Grafico de decis", xlab="Decis", ylab="Valores", ylim=c(mi,ma),
		type="b", col="red"); lines(q, q_y, type="b", col="blue")
	legend(0, ma, c("Grafico de decis 1", "Grafico de decis 2"), col=c("blue","red"), pch=rep(20,2))
	print(paste("Valores dos decis"))
	print(val)
}

#---------------------------------------------------------------------------
#TESTE UNIVARIADO PARA COMPARAcaO DE UM QUANTIL ESPECiFICO
#---------------------------------------------------------------------------
quniv.test.mqd<-function(x, med=med, q=.5, tipo="bilateral") {
	#Testa se o quantil q e igual a um valor hipotetico med.
	#O default e o teste bilateral para a mediana.
	x<-as.matrix(x)
	x<-exc.mv.mqd(x)
	x<-sort(x)
	n=length(x)
	amostral=round(hdquantil.mqd(x, q), 3)
	pval.inf=sum(dbinom(0:length(x[x<=med]), n, q))
	pval.sup=sum(dbinom((length(x[x<=med])+1):n, n, q))
	pval.bi=2*min(pval.inf, pval.sup)
	ic.med<-ic.q.mqd(x, q)
		
	if(tipo=="bilateral") {
		p_valor=round(pval.bi,3)
		ic=ic.med[2:3]
	}
	if(tipo=="menor") {
		p_valor=round(pval.sup,3)
		pos.sup=qbinom(.95, n, q)
		ls=round(x[pos.sup], 3)
		ic=c(LI="-inf", LS=ls)
	}
	if(tipo=="maior") {
		p_valor=round(pval.inf,3)
		pos.inf=qbinom(.05, n, q)
		li=round(x[pos.inf], 3)
		ic=c(LI=li, LS="+inf")
	}
	list(Q_suposto=med, Q_observado=amostral, IC=ic, p_valor=p_valor)
}
#---------------------------------------------------------------------------
#TESTE PARA POSIcaO CENTRAL
#---------------------------------------------------------------------------
mediana.mqd<-function(x,y, par="FALSE") { 
        if(par) {
	        dif=x-y
	        teste=quniv.test.mqd(dif, med=0)
	        res=list(p=teste$p_valor)
        }
        else {
	        z<-c(x,y) 
	        g <- rep(1:2, c(length(x),length(y)))
	        m<-median(z) 
	        teste=chisq.test(z<m,g)
		res=list(p=round(teste$p.value, 3)) 
        }
	res
}
#--------------------------------------------------------------------------
central.test.mqd<-function(x, y=NULL, par=FALSE, medida=0, tipo="bilateral") {
#Teste para a posicao central pelo teste t, teste de Wilcoxon e teste da mediana
#O default e o teste bilateral
	n=length(x)
	media=round(mean(x), 3)
	mediana=round(hdquantil.mqd(x, .5), 3)
	if(!is.null(y)) {
		n1=length(y)
		media_y=round(mean(y), 3)
		mediana_y=round(hdquantil.mqd(y, .5), 3)
	  	if(par==FALSE) {
   	  		print("Teste para duas amostras NaO pareadas")
	    		tt=t.test(x, y)
          		testemd=mediana.mqd(x, y, par=FALSE)
          		md.pval=testemd$p

		        testew=wilcox.test(x, y, paired=FALSE)
		        w.pval=round(testew$p.value, 3)
		        w.est=round(testew$statistic, 3) 
	  	}
		 if(par==TRUE) {
			print("Teste para duas amostras pareadas")
			tt=t.test(x, y, paired=TRUE)
		            
		        testemd=mediana.mqd(x, y, par=TRUE)
		        md.pval=testemd$p

		        testew=wilcox.test(x, y, paired=TRUE)
		        w.pval=round(testew$p.value, 3)
		        w.est=round(testew$statistic, 3) 

		}
	        print("Hipotese nula: a media e a mediana sao iguais")
	        print(paste("Teste t: Media1=", media, ", Media2=", media_y, ", Estatistica t=", round(tt$statistic,3), ",", n1+n-2, "gl", ", p=", round(tt$p.value,3)))   
	        print(paste("Teste da mediana: Med1=", mediana, ", Med2=", mediana_y, ", p_valor=", md.pval))
	        print(paste("Teste de Wilcoxon: Estatistica W=", w.est, ", p=", w.pval)) 
	}

 	else {
        	if(tipo=="bilateral") {
		        testet=t.test(x, mu=medida)
		        t.pval<-round(testet$p.value, 3)
		        t.est=round(testet$statistic, 3) 
		        li=testet$conf.int[1]; ls=testet$conf.int[2]
	          
		        testemd=quniv.test.mqd(x, medida, q=.5)
		        md.pval<-testemd$p_valor
		        md.para=testemd$IC
		        
		        testew=wilcox.test(x, mu=medida)
		        w.pval<-round(testew$p.value, 3)
		        w.est=round(testew$statistic, 3) 
	
		        print(paste("Hipotese nula: a media e a mediana sao iguais a", medida))
		        print(paste("Teste t: Media=", media, ", LI=", round(li,3), ", LS=", round(ls,3), "e p=", t.pval, ", com t=", t.est, "e", n-1, "gl"))   
		        print(paste("Teste da mediana: Med=", mediana, ", LI=", round(md.para[1],3), ", LS=", round(md.para[2],3), "e p=", md.pval))
		        print(paste("Teste de Wilcoxon: Med=", mediana, ", LI=", round(md.para[1],3), ", LS=", round(md.para[2],3), "e p=", w.pval, ", com V=", w.est))
        	}
        
	        if(tipo=="menor") {
		        testet=t.test(x, mu=medida, alternative="less")
		        t.pval<-round(testet$p.value, 3)
		        t.est=round(testet$statistic, 3) 
		        li=testet$conf.int[1]; ls=testet$conf.int[2]
		
		        testemd=quniv.test.mqd(x, medida, q=.5, "menor")
		        md.pval<-testemd$p_valor
		        md.para=testemd$IC
		
		        testew=wilcox.test(x, mu=medida, alternative="less")
		        w.pval<-round(testew$p.value, 3)
		        w.est=round(testew$statistic, 3) 
	
		        print(paste("Hipotese nula: a media e a mediana sao menores que", medida))
		        print(paste("Teste t: Media=", media, ", LI=", round(li,3), ", LS=", round(ls,3), "e p=", t.pval, ", com t=", t.est, "e", n-1, "gl"))   
		        print(paste("Teste da mediana: Med=", mediana, ", LI=", md.para[1], ", LS=", md.para[2], "e p=", md.pval))
		        print(paste("Teste de Wilcoxon: Med=", mediana, ", LI=", md.para[1], ", LS=", md.para[2], "e p=", w.pval, ", com V=", w.est))
	        }

	        if(tipo=="maior") {
		        testet=t.test(x, mu=medida, alternative="greater")
		        t.pval<-round(testet$p.value, 3)
		        t.est=round(testet$statistic, 3)
		        li=testet$conf.int[1]; ls=testet$conf.int[2]
		         
		        testemd=quniv.test.mqd(x, medida, q=.5, "maior")
		        md.pval<-testemd$p_valor
		        md.para=testemd$IC
		
		        testew=wilcox.test(x, mu=medida, alternative="greater")
		        w.pval<-round(testew$p.value, 3)
		        w.est=round(testew$statistic, 3) 
		        print(paste("Hipotese nula: a media e a mediana sao maiores que", medida))
		        print(paste("Teste t: Media=", media, ", LI=", round(li,3), ", LS=", round(ls,3), "e p=", t.pval, ", com t=", t.est, "e", n-1, "gl"))   
		        print(paste("Teste da mediana: Med=", mediana, ", LI=", md.para[1], ", LS=", md.para[2], "e p=", md.pval))
		        print(paste("Teste de Wilcoxon: Med=", mediana, ", LI=", md.para[1], ", LS=", md.para[2], "e p=", w.pval, ", com V=", w.est))
	        }
	}	
}

#---------------------------------------------------------------------------
#TESTE - TESTE PARA MEDIANA DE DUAS AMOSTRAS
#---------------------------------------------------------------------------
mediana.test.mqd<-function(x,y, par="FALSE") { 
	if(par) {
		print("Teste para duas amostras pareadas")
		dif=x-y
		teste=quniv.test.mqd(dif, med=0)
		res=c(Mediana_=c(hdquantil.mqd(x), hdquantil.mqd(y)), P_valor=teste$p_valor)
	}
	else {
		print("Teste para duas amostras independentes (com base no teste qui-quadrado)")
		z<-c(x,y) 
		g <- rep(1:2, c(length(x),length(y)))
		m<-median(z) 
		teste=chisq.test(z<m,g)
		x2=teste$statistic
		df=format(teste$parameter, scientific = F)
		res=c(Mediana_X=hdquantil.mqd(x), Mediana_Y=hdquantil.mqd(y), x2, df, P_valor=round(teste$p.value, 5))
	}
	res
} 

#---------------------------------------------------------------------------
#TESTE - PARTE 1: BOOSTRAPPING PARA UM QUANTIL ESPECiFICO
#---------------------------------------------------------------------------
qbiv.test.mqd=function(x, y, q=.5) {
	k=1000			#Numero de iteracoes
	t1=length(x)		#Tamanho da amostra 1
	t2=length(y)		#Tamanho da amostra 2
	Quantis=matrix(0, k, 3)
	for (i in 1:k) {
		Quantis[i,1]<-hdquantil.mqd(sample(x, t1, replace=TRUE), q)
		Quantis[i,2]<-hdquantil.mqd(sample(y, t2, replace=TRUE), q)
		Quantis[i,3]<-Quantis[i,1]-Quantis[i,2]
	}

	#Intervalo de confianca da diferenca de quantis
	q_x=round(hdquantil.mqd(x, q), 3); q_y=round(hdquantil.mqd(y, q), 3); q_d=round(q_x-q_y, 3)
	LI=round(hdquantil.mqd(Quantis[,3], 0.025), 3); LS=round(hdquantil.mqd(Quantis[,3], 0.975), 3)
	
	#Teste e hipoteses
	media=mean(Quantis[,3]); desvio=sd(Quantis[,3]); n=k; tcalc=(media)/(desvio)
	p_valor_1=round(2*min(pt(tcalc, n-1), 1-pt(tcalc, n-1)), 3)
	Resultado=ifelse(p_valor_1<0.05, "Rejeita a igualdade, a 5%", "Nao rejeitamos a igualdade, a 5%")

	p=(length(Quantis[,3][Quantis[,3]<0])+0.5*length(Quantis[,3][Quantis[,3]==0]))/k
	p_valor_2=round(2*min(p, 1-p), 3)

	print("--------------------------------------------------------------------")
	print("RESULTADOS")
	print(paste("O percentil do teste e:",100*q, "(o teste default e da mediana)"))
	print(paste("O quantil na amostra a e: ", q_x))
	print(paste("O quantil na amostra b e: ", q_y))
	print(paste("A diferenca de quantis e: ", q_d))
	print(paste("O limite inferior do IC de 95% e:", LI))
	print(paste("O limite superior do IC de 95% e:", LS))
	print(paste("O p-valor (calculado de suas formas) foi e: ", p_valor_1, "ou", p_valor_2))
	print(paste("O resultado do teste e:   ", Resultado))
	print("--------------------------------------------------------------------")
	print("DETALHES TeCNICOS")
	print("Estimadores de medidas por bootstrapping com 1000 iteracoes")
	print(paste("Hipotese nula: os quantis de ordem", q, "sao iguais"))
	print("Foram usados os estimadores de quantis de Harrell-Davis")
}

#---------------------------------------------------------------------------
#TESTE - PARTE 2: BOOSTRAPPING PARA DECIS
#---------------------------------------------------------------------------
decis.test.mqd=function(x, y) {
	k=1000			#Numero de iteracoes
	t1=length(x)		#Tamanho da amostra x
	t2=length(y)		#Tamanho da amostra y

	Dec=function(q) {
		Quantis=matrix(0, k, 3)
		for (i in 1:k){
			Quantis[i,1]<-hdquantil.mqd(sample(x, t1, replace=TRUE), q)
			Quantis[i,2]<-hdquantil.mqd(sample(y, t2, replace=TRUE), q)
			Quantis[i,3]<-Quantis[i,1]-Quantis[i,2]
		}
		LI=round(hdquantil.mqd(Quantis[,3], 0.025), 3); LS=round(hdquantil.mqd(Quantis[,3], 0.975), 3)
		p=(length(Quantis[,3][Quantis[,3]<0])+0.5*length(Quantis[,3][Quantis[,3]==0]))/k
		p_valor=round(2*min(p, 1-p), 3)
		c(LI, LS, p_valor)
	}

	Pos=matrix(0, 9, 7)
	for (i in 1:9) {
		Pos[i, 1]=i/10
		Pos[i, 2]=round(hdquantil.mqd(x, i/10), 3)
		Pos[i, 3]=round(hdquantil.mqd(y, i/10), 3)
		Pos[i, 4]=round(Pos[i, 2]-Pos[i, 3], 3)
		Pos[i, 5:7]=Dec(i/10)
	}
	plot(Pos[,1], Pos[,4], main="Grafico de decis", xlab="Decis", ylab="Valores",
	ylim=c(min(Pos[,5]), max(Pos[,6])), type="b"); 
	lines(Pos[,1], Pos[,5], type="b", col="blue");
	lines(Pos[,1], Pos[,6], type="b", col="blue");
	lines(Pos[,1], rep(0, 9),type="b", col="red");

	Vis=as.data.frame(Pos)
	names(Vis)=c("Decis", "Am 1", "Am 2", "Dif", "LI", "LS", "p-valor") 
	print("--------------------------------------------------------------------")
	print(Vis)
	print("--------------------------------------------------------------------")
	print("LEGENDAS")
	print("Decis: os 9 decis ordenados")
	print("Amos. 1 e 2: amostras 1 e 2")
	print("Dif.: diferenca dos quantis das duas amostras")
	print("LI e LS: limiites inferior e superior do intervalo de confianca de 95% para a diferenca")
	print("p-valor: nivel de significância para analise da hipotese nula")
	print("--------------------------------------------------------------------")
	print("DETALHES TeCNICOS")
	print("Estimadores de medidas por bootstrapping com 1000 iteracoes")
	print("A hipotese nula e de que nao ha diferenca entre os respectivos decis")
	print("Foram usados os estimadores de quantis de Harrell-Davis")
	print("O grafico: preto- diferenca; azul- limites do IC; vermelho - zero")

}

#---------------------------------------------------------------------------
#PARTE - TESTE GERAL PARA QUANTIS (AGREGA quniv e qbiv)
#---------------------------------------------------------------------------
q.test.mqd<-function(x, y=NULL, med=0, q=.5, tipo="bilateral") {
	if(is.null(y)) {
		quniv.test.mqd(x, med=med, q=.5, tipo="bilateral")
	}
	else {
		qbiv.test.mqd(x, y, q=q)
	}
}

#--------------------------------------------------------------
#TESTE DE CRAMER-VON MISES
#--------------------------------------------------------------
cvm.mqd<-function(x, y, plotit=TRUE) {
	N1=length(x); N2=length(y)
	N1=as.numeric(N1)
	N2=as.numeric(N2)
	N=N1+N2
	z=rank(c(x, y))             #Sobreponto
	aux=rep(1:2, c(N1, N2))     #Auxiliar
	cb=data.frame(c(x, y), aux,z)
	cb=cb[order(cb$z),]
	mx=data.frame(cb[which(cb$aux==1),], sq=seq(1:N1))
	my=data.frame(cb[which(cb$aux==2),], sq=seq(1:N2))
	sx=(mx$sq-mx$z)
	sy=(my$sq-my$z)
	U=N1*sum(sx^2)+N2*sum(sy^2)
	T=round((U/(N1*N2)/N)-(4*N1*N2-1)/(6*N), 3)

	#Referenciais para calcular o p-valor (baseado em CONOVER...)
	Tab=c(.046, .062, .079, .097, .119, .147, .184, .241, .347, .461, .743, 1.168)
	level=c(.1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .99, .999)
	
	#Graficos de funcoes de distribuicao
	if(plotit) {
		plot(ecdf(x)); lines(ecdf(y), col="red")
	}
	#O teste
	if(T>max(Tab)) p=0.000
	else if(T<min(Tab)) p=1.000
	else p=round(1-approx(Tab, level, T)$y, 3)
	
	Res <- list(Estatistica_T = T, p_valor = p)
	return(Res)
}

#===========================================================================
#     COMPARAcaO DE IGUALDADE DE DUAS DISTRIBUIcoES
#===========================================================================
dist.test.mqd<-function(x, y) {
	KS<-ks.test(x,y)
	CVM<-cvm.mqd(x,y, plotit=FALSE)
	print("Hipotese nula: as amostras sao oriundas de populacoes com distribuicoes idênticas")
	print(paste("Teste de Kolmogorov-Smirnov: Estatistica D=", round(KS$statistic, 3), ", p-valor=", round(KS$p.value, 3)))
	print(paste("Teste de Cramer-von Mises: Estatistica T=", CVM$Estatistica_T, ", p-valor=", CVM$p_valor))
}

#===========================================================================
#     COMPARAcaO DE QUANTIS DE AMOSTRAS INDEPENDENTES
#===========================================================================
#---------------------------------------------------------------------------
#PARTE - GERADOR DE MEDIDAS
#---------------------------------------------------------------------------
gera.quantis.mqd=function(x,y, quantidade=11) {
	y=as.factor(y)
	dados=data.frame(x,y)
	l=levels(y)
	g=length(l)
	Quantis=(0:(quantidade-1))/(quantidade-1)
	M=matrix(0, g, quantidade)
	for (i in 1:g) {
		M[i, 1:quantidade]=quantile(x[y==l[i]], Quantis)
		mi=min(M)	
		ma=max(M)
	}

	#Visualizacao grafica
	plot(Quantis, M[1,], main="Grafico de quantis", xlab="Quantis", ylab="Valores", 
			ylim=c(mi,ma), type="b", col=1); 

	for (i in 2:g) {
		lines(Quantis, M[i,], type="b", col=i)
	}

	print("--------------------------------------------------------------------")
	print("Funcao: gera.quantis.mqd(x=var. quantit., y=var. categ., quantidade=num. de quantis)")
	print(paste("Grafico de linhas e", quantidade, "quantis (default: minimo, maximo e 9 decis)"))
	print(paste("Os quantis dos", g, "grupos sao:"))
	round(data.frame(Quantis, t(M)),3)
}

#======================================================================================
#    CONTEuDOS DIVERSOS PARA MODELAGEM
#======================================================================================

#--------------------------------------------------------------------------------------
#DIAGPLOT PARA RESiDUOS
#--------------------------------------------------------------------------------------
diagplot.mqd<-function(modelo, com_out=TRUE, preditores=cbind(x), modelagem="lm") {
  	m = model.frame(modelo)[1]
	y = m[,1]
  	#Definindo os residuos
	if(modelagem %in% c("lm", "rfit")) {
	  	y_est=fitted(modelo)
	  	residuo = residuals(modelo)
  	}

	if(modelagem %in% c("tsreg","tshdreg", "tsregF", "tsregNW")) {
		y_est=est.mqd(modelo, preditores=preditores, modelagem=modelagem)
		residuo = y-y_est
	}

	respadron = (residuo-mean(residuo))/sd(residuo)
 
 	#Eliminando os outliers
 	y_sem_out<-exc.out.mqd(y)
	y_est_sem_out<-exc.out.mqd(y_est)
	sem_out<-data.frame(y_sem_out,y_est_sem_out)
	sem_out<-exc.mv.mqd(sem_out)
	
	#Definindo os residuos padronizados sem outliers
	res_sem_out<-sem_out[,1]-sem_out[,2]
	respadron_sem_out<-(res_sem_out-mean(res_sem_out))/sd(res_sem_out)

	if(com_out) {
		Res.padr=respadron
		Estimados=y_est
	}
	else {
		Res.padr=respadron_sem_out
		Estimados=sem_out[,1]
	}
	par(mfrow=c(2,2))
	plot(Estimados, Res.padr, main = "Residuos por Estimados"); abline(h = c(-2.5, 2.5))
	boxplot(Res.padr, main = "Box-plot dos residuos padronizados")
	par(mfrow=c(2,2))
	plot(Estimados, Res.padr, main = "Residuos por Estimados"); abline(h = c(-2.5, 2.5))
	boxplot(Res.padr, main = "Box-plot dos residuos padron.")
	hist(Res.padr, main = "Histograma dos residuos")
	env.mqd(respadron, Exibe=FALSE)
}


#--------------------------------------------------------------------------------------
#GERADOR DO ENVELOPE SIMULADO
#--------------------------------------------------------------------------------------
env.mqd=function(x, Exibe=TRUE, b=500) {
	n=length(x)
	x=sort(x)
	teo_prob=(seq_along(x)-.5)/n
	teo_quant=qnorm(teo_prob); y=teo_quant
	if(b==0) {
		quant=matrix(0, n, 3)
        	for (i in 1:n) {
        		quant[i, 1:3]=ic.q.mqd(teo_quant, i/n)
        	}
	}
	else {
		aux=replicate(b, sort(rnorm(n)))
		LI1=c(apply(aux, 1, quantile, probs=.025))
		LS1=c(apply(aux, 1, quantile, probs=.975))
		quant=round(cbind(Quantil=teo_quant, LI=LI1, LS=LS1), 3)
	}
	
	Resq=as.data.frame(round(quant, 3))
	names(Resq)=c("Valor", "LI", "LS") 
	
	mat=cbind(Resq[,2], Resq[,3])
	mat=as.matrix(mat)
	plot(teo_quant, x, main="Envelope normal", xlim=c(1.1*min(teo_quant), 1.1*max(teo_quant)), ylim=c(1.1*min(x), 1.1*max(x)))
	abline(lm(y~teo_quant))
	matlines(teo_quant, mat, lty=2, col="red")
	if(Exibe) Resq
}


#--------------------------------------------------------------------------------------
#TESTE DE JARQUE-BERA
#--------------------------------------------------------------------------------------

jarque.bera.mqd<-function(x) {
	n=length(x)
	s=(sum((x-mean(x))^3)/n)/sd(x)^3
	c=(sum((x-mean(x))^4)/n)/sd(x)^4
	est=round((n/6)*(s^2+0.25*(c-3)^2),3)
	pval<-round(1-pchisq(est, 2),3)
	list(Qui_quadrado=est, P_valor=pval)	
}


#--------------------------------------------------------------------------------------
#TESTES DE NORMALIDADE
#--------------------------------------------------------------------------------------
normalidade.mqd<-function(x) {
	cat("- Agrega seis testes de normalidade\n")
	cat("- Em todos, a hipotese nula e de que a distribuicao e normal\n")
	cat("\n")
	require(nortest)
	anderson<-ad.test(x)
	anderson.est<-round(anderson$statistic, 3)
	anderson.pval<-round(anderson$p.value, 3)
	anderson.res<-c("A", anderson.est, anderson.pval)
	
	lillie<-lillie.test(x)
	lillie.est<-round(lillie$statistic, 3)
	lillie.pval<-round(lillie$p.value, 3)
	lillie.res<-c("D", lillie.est, lillie.pval)
	
	shapiro<-sf.test(x)
	shapiro.est<-round(shapiro$statistic, 3)
	shapiro.pval<-round(shapiro$p.value, 3)
	shapiro.res<-c("W", shapiro.est, shapiro.pval)
	
	pearson<-pearson.test(x)
	pearson.est<-round(pearson$statistic, 3)
	pearson.pval<-round(pearson$p.value, 3)
	pearson.res<-c("P", pearson.est, pearson.pval)
	
	cramer<-cvm.test(x)
	cramer.est<-round(cramer$statistic, 3)
	cramer.pval<-round(cramer$p.value, 3)
	cramer.res<-c("W", cramer.est, cramer.pval)
	
	jarque<-jarque.bera.mqd(x)
	jarque.est<-jarque$Qui_quadrado
	jarque.pval<-jarque$P_valor
	jarque.res<-c("X²", jarque.est, jarque.pval)

	Resultado<-cbind(anderson.res, lillie.res, shapiro.res, pearson.res, cramer.res, jarque.res)
	Resultado=as.data.frame(Resultado)
	names(Resultado)=c("Anderson-Darling", "Lilliefors", "Shapiro-Francia", "Pearson X²",
			"Cramer-von Mises", "Jarque-Bera")
	res=as.data.frame(t(Resultado))
	names(res)=c("Estatistica", "Valor", "P.valor")
	res
}

#--------------------------------------------------------------------------------------
#TESTES DE HOMOSCEDASTICIDADE
#--------------------------------------------------------------------------------------
homoscedasticidade.mqd<-function(modelo, ncvTest=TRUE) {
	cat("- Agrega testes de homoscedasticidade\n")
	cat("- Em todos, a hipotese nula e de que os erros sao homoscedasticos\n")
	cat("\n")
	require(lmtest)
	require(car)
	BP<-bptest(modelo)
	BP.est<-round(BP$statistic, 3)
	BP.pval<-round(BP$p.value, 3)
	BP.res<-c("BP", BP.est, BP.pval)
	
	GQ<-gqtest(modelo)
	Est=as.numeric(format(GQ$statistic, scientific=F))
	GQ.est<-round(Est, 3)
	GQ.pval<-round(GQ$p.value, 3)
	GQ.res<-c("GQ", GQ.est, GQ.pval)

	if(ncvTest) {
		CW<-ncvTest(modelo)
		CW.est<-round(CW$ChiSquare, 3)
		CW.pval<-round(CW$p, 3)
		CW.res<-c("CW", CW.est, CW.pval)

		Resultado<-cbind(BP.res, GQ.res, CW.res)
		Resultado=as.data.frame(Resultado)
		names(Resultado)=c("Breush-Pagan", "Goldfeld-Quandt", "Cook-Weisberg")
 	}
	else {
		Resultado<-cbind(BP.res, GQ.res)
		Resultado=as.data.frame(Resultado)
		names(Resultado)=c("Breush-Pagan", "Goldfeld-Quandt")
 	}

	res=as.data.frame(t(Resultado))
	names(res)=c("Estatistica", "Valor", "P.valor")
	res
}

#--------------------------------------------------------------------------------------
#TESTES DE INDEPENDÊNCIA DE ERROS
#--------------------------------------------------------------------------------------
independencia.mqd<-function(modelo) {
	cat("- Agrega testes de independência de erros\n")
	cat("- Em todos, a hipotese nula e de que os erros sao independentes\n")
	cat("\n")
	require(lmtest)
	require(lawstat)
	
	DW<-dwtest(modelo)
	DW.est<-round(DW$statistic, 3)
	DW.pval<-round(DW$p.value, 3)
	DW.res<-c("DW", DW.est, DW.pval)
	
	BG<-bgtest(modelo)
	BG.est<-round(BG$statistic, 3)
	BG.pval<-round(BG$p.value, 3)
	BG.res<-c("BG", BG.est, BG.pval)
	
	RT<-runs.test(modelo$residuals)
	RT.est<-round(RT$statistic, 3)
	RT.pval<-round(RT$p.value, 3)
	RT.res<-c("Runs", RT.est, RT.pval)

	Resultado<-cbind(DW.res, BG.res, RT.res)
	Resultado=as.data.frame(Resultado)
	names(Resultado)=c("Durbin-Watson", "Breusch-Godfrey", "Teste de runs")
	res=as.data.frame(t(Resultado))
	names(res)=c("Estatistica", "Valor", "P.valor")
	res
}

#--------------------------------------------------------------------------------------
#MODELAGEM COMPLETA - NORMAL LINEAR
#--------------------------------------------------------------------------------------
mnl.mqd=function(formula, stepwise=FALSE) {
	require(MASS)
	modelo=lm(formula)
	sumario=summary(modelo)
	respad = ls.diag(modelo)$std.res
	diagplot.mqd(modelo)

	cat("MODELO, GRaFICOS E DIAGNOSTICO\n")
	cat("-----------------------MODELO-----------------------\n")
	coef=sumario$coefficients; ic=confint(modelo)
	result=round(cbind(coef[,1], coef[,1]-2*coef[,2], coef[,1]+2*coef[,2], coef[,2],coef[,3], coef[,4]),3)
	result=as.data.frame(result)
	names(result)=c("Estimad.", "LI (2.5%)", "LS (97.5%)", "Erro pad.", "Est. t", "p-valor")
	print(result)
	cat("\n")
	cat("--------------------AJUSTE GLOBAL--------------------\n")
	ag=sumario$fstatistic
	pval=1-pf(ag[1], ag[2], ag[3])
	cat(paste("Estatitica F=", round(ag[1], 3), ", gl=", ag[2], "e", ag[3], ", p-valor=", pval))
	cat("\n")
	r2=round(sumario$r.squared, 3); adj.r2=round(sumario$adj.r.squared, 3)
	cat(paste("R-quadrado=", r2, ", R-quadrado ajustado=", adj.r2, "\n"))
	cat("\n")

	cat("------------TESTES DE HOMOSCEDASTICIDADE------------\n")
	print(homoscedasticidade.mqd(modelo))
	cat("\n")
	
	cat("---------------TESTES DE NORMALIDADE----------------\n")
	print(normalidade.mqd(respad))
	cat("\n")
	
	cat("---------------TESTES DE INDEPENDÊNCIA--------------\n")
	print(independencia.mqd(modelo))
	cat("\n")

	if(stepwise) {
		cat("-----------MELHORIA DO MODELO - STEPWISE------------\n")
		step_modelo=stepAIC(modelo, direction="both")
		final=step_modelo$coefficients
		cat("\n")
		cat("- MODELO FINAL DO STEPWISE\n")
		print(final)
	}

}

#--------------------------------------------------------------------------------------
#INTERVALO DE CONFIANcA PARA O MODELO RFIT
#--------------------------------------------------------------------------------------
confint.rfit.mqd=function(modelo, nivel=.95) {
	alfa=1-nivel
	a=data.frame(summary(modelo)$coefficients)
	n_coef=length(a[,2])
	ic=matrix(0, n_coef, 3)
	for(i in 1:n_coef) {
		ic[i, 1:3]=round(c(a[i,1], a[i,1]+qnorm(alfa/2)*a[i,2], a[i,1]+qnorm(1-alfa/2)*a[i,2]), 3)
	}
	print(paste("Intervalo de confianca de", 100*nivel, "%"))
	print("LI e o limite inferior e LS e o limite superior do intervalo de confianca")
	ic=cbind(row.names(a),ic)
	ic=as.data.frame(ic)
	names(ic)=c("Variaveis", "Estimador", "LI", "LS")
	ic
}

#--------------------------------------------------------------------------------------
#MODELAGEM COMPLETA - MODELAGEM COM RANQUES
#--------------------------------------------------------------------------------------
mbr.mqd=function(formula) {
	require(Rfit)
	library(car)
	require(MASS)
	
	modelo=rfit(formula)
	sumario=summary(modelo)
	respad = rstudent(modelo)
	diagplot.mqd(modelo, modelagem="rfit")
	cat("MODELO, GRaFICOS E DIAGNOSTICO\n")
	cat("-----------------------MODELO-----------------------\n")
	coef=sumario$coefficients; ic=confint(modelo)
	coef=sumario$coefficients
	result=round(cbind(coef[,1], coef[,1]-2*coef[,2], coef[,1]+2*coef[,2], coef[,2],coef[,3], coef[,4]),3)
	result=as.data.frame(result)
	names(result)=c("Estimad.", "LI (2.5%)", "LS (97.5%)", "Erro pad.", "Est. t", "p-valor")
	print(result)

	cat("\n")
	cat("--------------------AJUSTE GLOBAL--------------------\n")
	ag=sumario$dropstat
	pval=sumario$droppval
	cat(paste("Teste de reducao da dispersao=", round(ag, 3), ", p-valor=", pval, "\n"))
	
	r2=round(sumario$R2, 3)
	cat(paste("R-quadrado robusto=", r2, "\n"))

	cat("\n")
	cat("------------TESTES DE HOMOSCEDASTICIDADE------------\n")
	print(homoscedasticidade.mqd(modelo, ncvTest=FALSE))
	cat("\n")
	
	cat("---------------TESTES DE NORMALIDADE----------------\n")
	print(normalidade.mqd(respad))
	cat("\n")
	
	cat("---------------TESTES DE INDEPENDÊNCIA--------------\n")
	print(independencia.mqd(modelo))

}

#--------------------------------------------------------------------------------------
#MODELAGEM COMPLETA - PSEUDO R2 QUANTiLICO
#--------------------------------------------------------------------------------------
mrqr2.mqd<-function(formula, quantil) {
	a=rq(formula, tau=quantil)				
	SPRA1=sum(abs(quantil*(y-fitted(a))[y>=fitted(a)]))
	SPRA2=sum(abs((1-quantil)*(y-fitted(a))[y<fitted(a)]))
	SPRA=SPRA1+SPRA2
	 
	SPTA1=sum(abs(quantil*(y-quantile(y,quantil))[y>=quantile(y,quantil)]))
	SPTA2=sum(abs((1-quantil)*(y-quantile(y,quantil))[y<quantile(y,quantil)]))
	SPTA=SPTA1+SPTA2
	FalsoR2=1-SPRA/SPTA
	cat("Pseudo_R2=")
	cat(round(FalsoR2, 3),"\n")
}

#--------------------------------------------------------------------------------------
#MODELAGEM COMPLETA - MODELO QUANTiLICO
#--------------------------------------------------------------------------------------
mrq_aux<-function(formula, q, se="boot") {
	require(quantreg)
	modelo=rq(formula, tau=q)
	sumario1=summary(modelo)
	sumario2=summary(modelo, se=se)

	cat(paste("------------MODELO PARA O QUANTIL", q, "-------------\n"))
	coef1=sumario1$coefficients
	coef2=sumario2$coefficients
	result=round(cbind(coef1, coef2[,2],coef2[,3], coef2[,4]),3)
	result=as.data.frame(result)
	names(result)=c("Estimad.", "LI (2.5%)", "LS (97.5%)", "Erro pad.", "Est. t", "p-valor")
	quantil=c("Quantil=", sumario1$tau)
	print(result)

	cat("\n")
	cat("--------------------AJUSTE GLOBAL--------------------\n")
	print(mrqr2.mqd(formula, q))
}

mrq.plot<-function(formula, q) {
	modelo=summary(rq(formula, tau=q))
	plot(modelo)
}

mrq.mqd<-function(formula, q=.5, se="boot") {
	print("Opcoes de erro padrao: boot (default), iid, nid")	
	for(i in 1:length(q)) {
		print(mrq_aux(formula, q[i], se=se))
	}
	cat("\n")
}

#======================================================================================
#    SIMULAcaO DE ESTIMADORES DE PARÂMETROS E INTERVALOS DE CONFIANcA
#======================================================================================
#--------------------------------------------------------------------------------------
#ANOVA ONEWAY
#--------------------------------------------------------------------------------------
oneway.mqd<-function(resp, categ) {
#Gera a anova parametrica e nao parametrica (Kruskal-Wallis)
#Ordenadamente, temos a variavel resposta e a categorica

	if(is.numeric(resp)==FALSE) stop("A variavel resposta precisa ser quantitativa")
	if(is.factor(categ)==FALSE) stop("A segunda variavel precisa ser categorica")
	aov.par<-summary(aov(resp~categ))
	est.par=round(aov.par[[1]]$'F value'[1], 3)
	df.par=aov.par[[1]]$Df
	pval.par=round(aov.par[[1]]$'Pr(>F)'[1], 3)
	
	aov.kw<-kruskal.test(resp~categ)
	est.kw=round(aov.kw$statistic, 3)
	df.kw=c(aov.kw$parameter, '-')
	pval.kw=round(aov.kw$p.value, 3)
	resultado=cbind(c(Estatistica=est.par, gl=df.par, P_valor=pval.par),
	c(Estatistica=est.kw, gl=df.kw, P_valor=pval.kw))
	resultado<-data.frame(resultado)
	names(resultado)=c("Anova classica", "Anova KW")
	resultado
}

#--------------------------------------------------------------------------------------
#MEDIDAS DESCRITIVAS DE POSIcaO
#--------------------------------------------------------------------------------------
descritivo.oneway.mqd<-function(resp, categ) {
	#Gera as medidas descritivas da variavel resposta por categoria
 
	aov.par=round(tapply(resp, categ, mean), 3)
	aov.kw1=round(tapply(resp, categ, median), 3)
	aov.kw2=round(tapply(rank(resp), categ, mean), 3)

	resultado=cbind(aov.par, aov.kw1, aov.kw2)
	resultado<-data.frame(resultado)
	colnames(resultado)="Grupos"
	names(resultado)=c("Medias", "Medianas", "Medias_postos")
	print(resultado)
}

#--------------------------------------------------------------------------------------
#TESTES DE COMPARcaO DOIS A DOIS
#--------------------------------------------------------------------------------------
posthoc.oneway.mqd<-function(resp, categ) {
#Gera as comparacoes dois a dois, para a anova classica e nao parametrica (Kruskal-Wallis)

	require(utils)
	cats <- split(resp, categ)
	nomes <- names(cats)
	combina=combn(nomes, 2)
	comps=ncol(combina)
	resultado<-matrix(0, comps, 3)
	for(i in 1:comps) {
		i1=combina[1,i]
		i2=combina[2,i]
		teste=wilcox.test(resp[categ==i1 | categ==i2]~categ[categ==i1 | categ==i2])
		est=teste$statistic
		pval=teste$p.value
		resultado[i, 1]=i1
		resultado[i, 2]=i2
		resultado[i, 3]=round(pval, 3)
	}
	resultado<-data.frame(resultado)
	names(resultado)=c("Grupo 1", "Grupo 2", "P_Wilxoxon")
	tuk=TukeyHSD(aov(resp~categ))$categ
	res=cbind(resultado, P_Tukey=c(round(tuk[,4], 3)))
	as.matrix(res)
}

#--------------------------------------------------------------------------------------
#RESULTADOS GERAS
#--------------------------------------------------------------------------------------
geral.oneway.mqd<-function(resp, categ) {
	require(Rfit)
	print("ANOVA CLaSSICA E DE KRUSKAL-WALLIS")
	print("___________________________________________________________")
	
	print(".........................TESTE.............................")
	print(oneway.mqd(resp, categ))
	
	print("..................COMPARAcaO DOIS A DOIS...................")
	print(posthoc.oneway.mqd(resp, categ))
	
	print("....................MEDIDAS DESCRITIVAS....................")
	print(descritivo.oneway.mqd(resp, categ))

	print("___________________________________________________________")
	print(".....................ANOVA POR RANQUES.....................")
	print(oneway.rfit(resp, categ))
	
	print("..................COMPARAcaO DOIS A DOIS...................")
	print(summary(oneway.rfit(resp, categ)))
	print("___________________________________________________________")
}

#======================================================================================
#    SIMULAcaO DE ESTIMADORES DE PARÂMETROS E INTERVALOS DE CONFIANcA
#======================================================================================
#----------------------------------------------------------------------------
#TESTE DE INTERVALOS DE CONFIANcA PARA A DISTRIBUIcaO DE BERNOULLI
#----------------------------------------------------------------------------
#Funcao: simula.bern.mqd(p, n, b, percentual5)
	#p - Parâmetro da distribuicao de Bernoulli
	#b - quantidade de reamostragens  (default de 3000) 
	#n - Tamanho da amostra (default de 50)
	#percentual - tamanho do intervalo de confianca (default de 95%)
#----------------------------------------------------------------------------
simula.bern.mqd = function(p, n=50, b=3000, percentual=95) {
	print("*********************************************")
	print(paste("MOSTRA INTERVALOS PARA O PARÂMETRO p DE X, X~Bernoulli(p)"))
	print("Funcao (defualt): IC.unif.mqd(p, tam. amostra=50, amostras=3000, % do IC=95)")	
	print(paste("Parâmetro da distribuicao:", p))
	print(paste("Tamanho de cada amostra (default de 50):", n))
	print(paste("Quantidade de amostras (default):", b))
	print(paste("Tamanho do intervalo (default):", percentual,"%"))	
	alfa=(100-percentual)/100
	z=qnorm(1-alfa/2)	

	M = matrix(0, b, n+14)
	for (i in 1:b) {
		set.seed(i)
		#Primeiro intervalo - classico
		M[i,1:n]=rbinom(n, 1, p)
		M[i,n+1]=mean(M[i,1:n])+qt(alfa/2, n-1)*((mean(M[i,1:n])*(1-mean(M[i,1:n])))/n)^0.5
		M[i,n+2]=mean(M[i,1:n])+qt(1-alfa/2, n-1)*((mean(M[i,1:n])*(1-mean(M[i,1:n])))/n)^0.5			
		M[i,n+3]=M[i,n+1]<p & M[i,n+2]>p

		#Segundo intervalo - Wilson
		M[i,n+4]=(1/(1+z^2/n))*(mean(M[i,1:n])+(z^2)/(2*n) - z*(((mean(M[i,1:n])*(1-mean(M[i,1:n])))/n)+(z^2)/(4*n^2))^0.5)
		M[i,n+5]=(1/(1+z^2/n))*(mean(M[i,1:n])+(z^2)/(2*n) + z*(((mean(M[i,1:n])*(1-mean(M[i,1:n])))/n)+(z^2)/(4*n^2))^0.5)
		M[i,n+6]=M[i,n+4]<p & M[i,n+5]>p

		#Terceiro intervalo - Agresti-Coull
		n_til=n+z^2
		p_til=(sum(M[i,1:n])+z^2/2)/n_til
		M[i,n+7]=p_til-z*(p_til*(1-p_til)/n)^0.5
		M[i,n+8]=p_til+z*(p_til*(1-p_til)/n)^0.5
		M[i,n+9]=M[i,n+7]<p & M[i,n+8]>p

		#Quarto intervalo - Arco_seno
		M[i,n+10]=(sin(asin(mean(M[i,1:n])^0.5)-z/(2*n^0.5)))^2
		M[i,n+11]=(sin(asin(mean(M[i,1:n])^0.5)+z/(2*n^0.5)))^2
		M[i,n+12]=M[i,n+10]<p & M[i,n+11]>p

	}

	Media_LI_1 = round(mean(M[,n+1]),3) 
	Media_LS_1 = round(mean(M[,n+2]),3) 
	Acerto_perc._1 = 100*round(mean(M[,n+3]),3)
	
	Media_LI_2 = round(mean(M[,n+4]),3) 
	Media_LS_2 = round(mean(M[,n+5]),3) 
	Acerto_perc._2 = 100*round(mean(M[,n+6]),3)

	Media_LI_3 = round(mean(M[,n+7]),3) 
	Media_LS_3 = round(mean(M[,n+8]),3) 
	Acerto_perc._3 = 100*round(mean(M[,n+9]),3)

	Media_LI_4 = round(mean(M[,n+10]),3) 
	Media_LS_4 = round(mean(M[,n+11]),3) 
	Acerto_perc._4 = 100*round(mean(M[,n+12]),3)

	print("-----------------------------------------------------")
	print(paste(" ** ANaLISE DO INTERVALO CLaSSICO ** "))
	print(paste("Media do limite inferior do IC: ", Media_LI_1))
	print(paste("Media do limite superior do IC: ", Media_LS_1))
	print(paste("Acerto percentual do IC: ", Acerto_perc._1, "%"))

	print("-----------------------------------------------------")
	print(paste(" ** ANaLISE DO INTERVALO DE WILSON** "))
	print(paste("Media do limite inferior do IC: ", Media_LI_2))
	print(paste("Media do limite superior do IC: ", Media_LS_2))
	print(paste("Acerto percentual do IC: ", Acerto_perc._2, "%"))

	print("-----------------------------------------------------")
	print(paste(" ** ANaLISE DO INTERVALO DE AGRESTI-COULL** "))
	print(paste("Media do limite inferior do IC: ", Media_LI_3))
	print(paste("Media do limite superior do IC: ", Media_LS_3))
	print(paste("Acerto percentual do IC: ", Acerto_perc._3, "%"))

	print("-----------------------------------------------------")
	print(paste(" ** ANaLISE DO INTERVALO ARCO-SENO** "))
	print(paste("Media do limite inferior do IC: ", Media_LI_4))
	print(paste("Media do limite superior do IC: ", Media_LS_4))
	print(paste("Acerto percentual do IC: ", Acerto_perc._4, "%"))

	print("-----------------------------------------------------")
	print("DETALHES")
	print(paste("IC classico: [media (- E +) t(alfa/2, n-1)*desvio/sqrt(n)]"))
	print(paste("IC de Wilson: [(1/(1+z^2/n)*media+(z^2/(2*n)) (- E +) z*(desvio/n+z^2/(4*sqrt(n))), com z=qnorm(1-alfa/2)]"))
	print(paste("IC de Agresti-Coull: [pt (- E +) z*(pt*(1-pt)/n)^0.5, com nt=n+z^2, pt=(soma+z^2/2)/nt, z=qnorm(1-alfa/2)]"))
	print(paste("IC de Arco_seno: [(sen(arcsen(sqrt(media)) (- E +) z/(2*sqrt(n))))^2, z=qnorm(1-alfa/2)]"))
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#Funcao: simula.unif.mqd(p, n, b, percentual5)
	#min - Parâmetro do minimo
	#max - Parâmetro do maximo
	#b - quantidade de reamostragens (default de 3000) 
	#n - Tamanho da amostra (default de 50)
	#percentual - tamanho do intervalo de confianca (default de 95%)

simula.unif.mqd=function(min, max,  n=50, b=3000, percentual=95) {
	print("*********************************************")
	print(paste("MOSTRA INTERVALOS PARA MiNIMO, MaXIMO E MeDIA DE X, X~U(min, max)"))
	print("Funcao (defualt): simula.unif.mqd(min, max, tam. amostra=50, amostras=3000, % do IC=95)")	
	print(paste("Parâmetros da distribuicao:", min, "e", max))
	print(paste("Valor esperado (media):", round((min+max)/2,2)))
	print(paste("Tamanho de cada amostra (default de 50):", n))
	print(paste("Quantidade de amostras (default):", b))
	print(paste("Tamanho do intervalo (default):", percentual,"%"))
	alfa=(100-percentual)/100
	gama=exp(-(log(alfa))/(n-1))-1
	media=(min+max)/2
	M = matrix(0, b, n+14)
	for (i in 1:b) {
		set.seed(i)
		M[i,1:n]=runif(n, min, max)
		M[i,n+1]=max-(max-min(M[i,1:n]))/(alfa/2)^(1/n)		#Supoe-se saber o maximo	
		M[i,n+2]=max-(max-min(M[i,1:n]))/(1-alfa/2)^(1/n)			
		M[i,n+3]=M[i,n+1]<min & M[i,n+2]>min

		M[i,n+4]=min+(max(M[i,1:n])-min)/(1-alfa/2)^(1/n)	#Supoe-se saber o mminimo	
		M[i,n+5]=min+(max(M[i,1:n])-min)/(alfa/2)^(1/n)		
		M[i,n+6]=M[i,n+4]<max & M[i,n+5]>max
	
		M[i,n+7]=(max(M[i,1:n])+min(M[i,1:n]))
		M[i,n+8]=(max(M[i,1:n])-min(M[i,1:n]))
		M[i,n+9]=0.5*(M[i,n+7]-M[i,n+8]*gama)
		M[i,n+10]=0.5*(M[i,n+7]+M[i,n+8]*gama)
		M[i,n+11]=M[i,n+9]<media & M[i,n+10]>media

		M[i,n+12]=mean(M[i,1:n])+qt(alfa/2, n-1)*sd(M[i,1:n])/n^0.5
		M[i,n+13]=mean(M[i,1:n])+qt(1-alfa/2, n-1)*sd(M[i,1:n])/n^0.5
		M[i,n+14]=M[i,n+12]<media & M[i,n+13]>media

	}
	
	Media_LI_Min = round(mean(M[,n+1]),3) 
	Media_LS_Min = round(mean(M[,n+2]),3) 
	Acerto_perc._Min = 100*round(mean(M[,n+3]),3)

	Media_LI_Max = round(mean(M[,n+4]),3) 
	Media_LS_Max = round(mean(M[,n+5]),3) 
	Acerto_perc._Max = 100*round(mean(M[,n+6]),3)
	
	Media_LI_Media_1 = round(mean(M[,n+9]),3) 
	Media_LS_Media_1 = round(mean(M[,n+10]),3) 
	Acerto_perc._Media_1 = 100*round(mean(M[,n+11]),3)

	Media_LI_Media_2 = round(mean(M[,n+12]),3) 
	Media_LS_Media_2 = round(mean(M[,n+13]),3) 
	Acerto_perc._Media_2 = 100*round(mean(M[,n+14]),3)
	
	print("-----------------------------------------------------")
	print(paste(" ** ANaLISE DO INTERVALO DE CONFIANcA DO MiNIMO ** "))
	print(paste("Media do limite inferior do IC do Minimo: ", Media_LI_Min))
	print(paste("Media do limite superior do IC do Minimo: ", Media_LS_Min))
	print(paste("Acerto percentual do IC: ", Acerto_perc._Min, "%"))

	print("-----------------------------------------------------")
	print(paste(" ** ANaLISE DO INTERVALO DE CONFIANcA DO MaXIMO ** "))
	print(paste("Media do limite inferior do IC do Maximo: ", Media_LI_Max))
	print(paste("Media do limite superior do IC do Maximo: ", Media_LS_Max))
	print(paste("Acerto percentual do IC: ", Acerto_perc._Max, "%"))

	print("-----------------------------------------------------")	
	print(paste(" ** ANaLISE DO PRIMEIRO INTERVALO DE CONFIANcA DA MeDIA ** "))
	print(paste("Media do limite inferior do IC_1 da Media: ", Media_LI_Media_1))
	print(paste("Media do limite superior do IC_1 da Media: ", Media_LS_Media_1))
	print(paste("Acerto percentual do IC: ", Acerto_perc._Media_1, "%"))

	print(paste(" ** ANaLISE DO SEGUNDO INTERVALO DE CONFIANcA DA MeDIA ** "))
	print(paste("Media do limite inferior do IC_2 da Media: ", Media_LI_Media_2))
	print(paste("Media do limite superior do IC_2 da Media: ", Media_LS_Media_2))
	print(paste("Acerto percentual do IC: ", Acerto_perc._Media_2, "%"))

	print("-----------------------------------------------------")	
	print("DETALHES")
	print(paste("IC do minimo: [max-(max-min(X))/(alfa/2)^(1/n);max-(max-min(X))/(1-alfa/2)^(1/n)], supoe-se saber o maximo"))
	print(paste("IC do maximo: [min+(max(X)-min)/(1-alfa/2)^(1/n);min+(max(X)-min)/(alfa/2)^(1/n)], supoe-se saber o minimo"))
	print(paste("Segundo IC da mediao: [0.5*((max(X)+min(X)) (- E +) (max(X)-min(X))*gama)], com gama=exp(-(log(alfa))/(n-1))-1"))
	print(paste("Segundo IC da mediao: [media (- E +) t(alfa/2, n-1)*desvio/sqrt(n)]"))
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
simula.pois.mqd = function(lambda, n=50, b=3000, percentual=95) {
	print("*********************************************")
	print(paste("MOSTRA INTERVALOS PARA O PARÂMETRO lambda DE X, X~Poisson(lambda)"))
	print("Funcao (defualt): simula.pois.mqd(lambda, tam. amostra=50, amostras=3000, % do IC=95)")	
	print(paste("Parâmetro da distribuicao:", lambda))
	print(paste("Tamanho de cada amostra (default de 50):", n))
	print(paste("Quantidade de amostras (default):", b))
	print(paste("Tamanho do intervalo (default):", percentual,"%"))	
	alfa=(100-percentual)/100

	M = matrix(0, b, n+11)
	for (i in 1:b) {
		set.seed(i)
		M[i,1:n] = rpois(n, lambda)
		M[i,n+1] = mean(M[i,1:n])
		M[i,n+2] = sd(M[i,1:n])
		M[i,n+3] = M[i,n+1]+qt(alfa/2,n-1)*M[i,n+2]/n^0.5
		M[i,n+4] = M[i,n+1]+qt(1-alfa/2, n-1)*M[i,n+2]/n^0.5
		M[i,n+5] = (lambda >= M[i,n+3]) & (lambda <= M[i,n+4])

		M[i,n+6] = (((qnorm(1-alfa/2)/n^0.5)-((qnorm(1-alfa/2)/n)+(4*M[i,n+1]))^0.5)/2)^2
		M[i,n+7] = (((qnorm(1-alfa/2)/n^0.5)+((qnorm(1-alfa/2)/n)+(4*M[i,n+1]))^0.5)/2)^2
		M[i,n+8] = (lambda >= M[i,n+6]) & (lambda <= M[i,n+7])

		M[i,n+9] = (qchisq(alfa/2, 2*n*M[i,n+1]))/(2*n)
		M[i,n+10] = (qchisq(1-alfa/2,2*(n*M[i,n+1]+2)))/(2*n)
		M[i,n+11] = (lambda >= M[i,n+9]) & (lambda <= M[i,n+10])
	}
	
	IC_Classico = 100*round(sum(M[,n+5])/b, 2)
	IC_Alternativo_1 = 100*round(sum(M[,n+8])/b, 2) 
	IC_Alternativo_2 = 100*round(sum(M[,n+11])/b, 2) 
	Media_Li_classico = round(mean(M[ ,n+3]),2)
	Media_Ls_classico = round(mean(M[ ,n+4]),2)
	Media_Li_alt_1 = round(mean(M[ ,n+6]),2)
	Media_Ls_alt_1 = round(mean(M[ ,n+7]),2)
	Media_Li_alt_2 = round(mean(M[ ,n+9]),2)
	Media_Ls_alt_2 = round(mean(M[ ,n+10]),2)
	
	print("-----------------------------------------------------")	
	print(paste(" ** ANaLISE DO PRIMEIRO INTERVALO DE CONFIANcA DA MeDIA (CLaSSICO) ** "))
	print(paste("Media do limite inferior do IC_2 da Media: ", Media_Li_classico))
	print(paste("Media do limite superior do IC_2 da Media: ", Media_Ls_classico))
	print(paste("Acerto percentual do IC: ", IC_Classico, "%"))
	print("-----------------------------------------------------")	
	print(paste(" ** ANaLISE DO SEGUNDO INTERVALO DE CONFIANcA DA MeDIA ** "))
	print(paste("Media do limite inferior do IC_2 da Media: ", Media_Li_alt_1))
	print(paste("Media do limite superior do IC_2 da Media: ", Media_Ls_alt_1))
	print(paste("Acerto percentual do IC: ", IC_Alternativo_1, "%"))
	print("-----------------------------------------------------")	
	print(paste(" ** ANaLISE DO TERCEIRO INTERVALO DE CONFIANcA DA MeDIA ** "))
	print(paste("Media do limite inferior do IC_2 da Media: ", Media_Li_alt_2))
	print(paste("Media do limite superior do IC_2 da Media: ", Media_Ls_alt_2))
	print(paste("Acerto percentual do IC: ", IC_Alternativo_2, "%"))
	print("-----------------------------------------------------")	
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
simula.norm.mqd = function(mu, sigma2, n=50, b=3000, percentual=95) {
	print("*********************************************")
	print(paste("MOSTRA INTERVALOS PARA OS PARÂMETROS mu E sigma² DE X, X~N(mu, sigma²)"))
	print("Funcao (defualt): simula.norm.mqd(mu, sigma2, tam. amostra=50, amostras=3000, % do IC=95)")	
	print(paste("Parâmetros da distribuicao: media -", mu, "variância -", sigma2))
	print(paste("Tamanho de cada amostra (default de 50):", n))
	print(paste("Quantidade de amostras (default):", b))
	print(paste("Tamanho do intervalo (default):", percentual,"%"))	
	alfa=(100-percentual)/100
	z=qnorm(1-alfa/2)
	ri=n/2-z*(n^0.5)/2
	rs=1+n/2+z*(n^0.5)/2
	
	M = matrix(0, b, n+14)
	for (i in 1:b) {
		set.seed(i)
		#Primeiro intervalo mu - classico
		M[i,1:n]=rnorm(n, mu, sigma2^0.5)
		M[i,n+1]=mean(M[i,1:n])-z*(sigma2/n)^0.5
		M[i,n+2]=mean(M[i,1:n])+z*(sigma2/n)^0.5
		M[i,n+3]=M[i,n+1]<mu & M[i,n+2]>mu
		
		#Segundo intervalo mu - baseado na distribuicao t		
		M[i,n+4]=mean(M[i,1:n])+qt(alfa/2, n-1)*(sd(M[i,1:n])/(n^0.5))
		M[i,n+5]=mean(M[i,1:n])+qt(1-alfa/2, n-1)*(sd(M[i,1:n])/(n^0.5))
		M[i,n+6]=M[i,n+4]<mu & M[i,n+5]>mu

		#Terceiro intervalo de mu - baseado na mediana
		M[i,n+7]=quantile(M[i,1:n], ri/n)
		M[i,n+8]=quantile(M[i,1:n], rs/n)
		M[i,n+9]=M[i,n+7]<mu & M[i,n+8]>mu

		#Terceiro intervalo de mu- baseado na mediana
		M[i,n+10]=((n-1)*var(M[i,1:n]))/qchisq(1-alfa/2, n-1) 
		M[i,n+11]=((n-1)*var(M[i,1:n]))/qchisq(alfa/2, n-1) 
		M[i,n+12]=M[i,n+10]<sigma2 & M[i,n+11]>sigma2
	}

	Media_LI_1 = round(mean(M[,n+1]),3) 
	Media_LS_1 = round(mean(M[,n+2]),3) 
	Acerto_perc._1 = 100*round(mean(M[,n+3]),3)

	Media_LI_2 = round(mean(M[,n+4]),3) 
	Media_LS_2 = round(mean(M[,n+5]),3) 
	Acerto_perc._2 = 100*round(mean(M[,n+6]),3)
	
	Media_LI_3 = round(mean(M[,n+7]),3) 
	Media_LS_3 = round(mean(M[,n+8]),3) 
	Acerto_perc._3 = 100*round(mean(M[,n+9]),3)

		
	Media_LI_4 = round(mean(M[,n+10]),3) 
	Media_LS_4 = round(mean(M[,n+11]),3) 
	Acerto_perc._4 = 100*round(mean(M[,n+12]),3)

	print("-----------------------------------------------------")
	print(paste(" ** ANaLISE DO INTERVALO DA MeDIA CLaSSICO (VARIÂNCIA CONHECIDA) ** "))
	print(paste("Media do limite inferior do IC: ", Media_LI_1))
	print(paste("Media do limite superior do IC: ", Media_LS_1))
	print(paste("Acerto percentual do IC: ", Acerto_perc._1, "%"))

	print(paste(" ** ANaLISE DO INTERVALO  DA MeDIA CLaSSICO (VARIÂNCIA DESCONHECIDA) ** "))
	print(paste("Media do limite inferior do IC: ", Media_LI_2))
	print(paste("Media do limite superior do IC: ", Media_LS_2))
	print(paste("Acerto percentual do IC: ", Acerto_perc._2, "%"))

	print(paste(" ** ANaLISE DO INTERVALO DA MeDIA BASEADO NA MEDIANA ** "))
	print(paste("Media do limite inferior do IC: ", Media_LI_3))
	print(paste("Media do limite superior do IC: ", Media_LS_3))
	print(paste("Acerto percentual do IC: ", Acerto_perc._3, "%"))

	print(paste(" ** ANaLISE DO INTERVALO DA VARIÂNCIA ** "))
	print(paste("Media do limite inferior do IC: ", Media_LI_4))
	print(paste("Media do limite superior do IC: ", Media_LS_4))
	print(paste("Acerto percentual do IC: ", Acerto_perc._4, "%"))
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
simula.exp.mqd = function(lambda, n=50, b=3000, percentual=95) {
	print("*********************************************")
	print(paste("MOSTRA INTERVALOS PARA A MeDIA DE X, X~Exp(lambda)"))
	print("Funcao (defualt): simula.mqd(lambda, tam. amostra=50, amostras=3000, % do IC=95)")	
	print(paste("Parâmetro da distribuicao:", round(lambda,2), "e a media e", round(1/lambda,2)))
	print(paste("Tamanho de cada amostra (default de 50):", n))
	print(paste("Quantidade de amostras (default):", b))
	print(paste("Tamanho do intervalo (default):", percentual,"%"))	
	alfa=(100-percentual)/100

	M = matrix(0, b, n+11)
	for (i in 1:b) {
		set.seed(i)
		M[i,1:n] = rexp(n, lambda)
		M[i,n+1] = mean(M[i,1:n])
		M[i,n+2] = sd(M[i,1:n])
		M[i,n+3] = M[i,n+1]+qt(alfa/2,n-1)*M[i,n+2]/n^0.5
		M[i,n+4] = M[i,n+1]+qt(1-alfa/2, n-1)*M[i,n+2]/n^0.5
		M[i,n+5] = (1/lambda >= M[i,n+3]) & (1/lambda <= M[i,n+4])

		M[i,n+6] = (n^.5*M[i,n+1])/(n^.5+qnorm(1-alfa/2))
		M[i,n+7] = (n^.5*M[i,n+1])/(n^.5-qnorm(1-alfa/2))
		M[i,n+8] = (1/lambda >= M[i,n+6]) & (1/lambda <= M[i,n+7])

		M[i,n+9] = (2*n*M[i,n+1])/qchisq(1-alfa/2,2*n)
		M[i,n+10] = (2*n*M[i,n+1])/qchisq(alfa/2,2*n)
		M[i,n+11] = (1/lambda >= M[i,n+9]) & (1/lambda <= M[i,n+10])
	}
	
	IC_Classico = 100*round(sum(M[,n+5])/b, 3)
	IC_Alternativo_1 = 100*round(sum(M[,n+8])/b, 3) 
	IC_Alternativo_2 = 100*round(sum(M[,n+11])/b, 3) 
	Media_Li_classico = round(mean(M[ ,n+3]),2)
	Media_Ls_classico = round(mean(M[ ,n+4]),2)
	Media_Li_alt_1 = round(mean(M[ ,n+6]),2)
	Media_Ls_alt_1 = round(mean(M[ ,n+7]),2)
	Media_Li_alt_2 = round(mean(M[ ,n+9]),2)
	Media_Ls_alt_2 = round(mean(M[ ,n+10]),2)
	
	print("-----------------------------------------------------")	
	print(paste(" ** ANaLISE DO PRIMEIRO INTERVALO DE CONFIANcA DA MeDIA (CLaSSICO) ** "))
	print(paste("Media do limite inferior do IC_2 da Media: ", Media_Li_classico))
	print(paste("Media do limite superior do IC_2 da Media: ", Media_Ls_classico))
	print(paste("Acerto percentual do IC: ", IC_Classico, "%"))
	print("-----------------------------------------------------")	
	print(paste(" ** ANaLISE DO SEGUNDO INTERVALO DE CONFIANcA DA MeDIA ** "))
	print(paste("Media do limite inferior do IC_2 da Media: ", Media_Li_alt_1))
	print(paste("Media do limite superior do IC_2 da Media: ", Media_Ls_alt_1))
	print(paste("Acerto percentual do IC: ", IC_Alternativo_1, "%"))
	print("-----------------------------------------------------")	
	print(paste(" ** ANaLISE DO SEGUNDO INTERVALO DE CONFIANcA DA MeDIA ** "))
	print(paste("Media do limite inferior do IC_2 da Media: ", Media_Li_alt_2))
	print(paste("Media do limite superior do IC_2 da Media: ", Media_Ls_alt_2))
	print(paste("Acerto percentual do IC: ", IC_Alternativo_2, "%"))
	print("-----------------------------------------------------")	
}

