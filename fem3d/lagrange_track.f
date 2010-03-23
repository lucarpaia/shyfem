c
c $Id: lagrange_track.f,v 1.1 2009-02-13 17:22:44 georg Exp $
c
c routines for tracking particle
c
c revision log :
c
c 05.02.2009    ggu     copied from other files
c
c**********************************************************************

        subroutine track_orig(time,bdy,ie,xbdy,ybdy)

c in questa subroutine si calcola il percorso del body partendo
c da un punto all'interno dell'elemento ie 
c fino al successivo punto interno allo stesso elemento
c se entro il time step il body esce dall'elemento 
c si passa alla successiva subroutine con time

        implicit none
        
	include 'param.h' 
	include 'lagrange.h'
 
	integer ie,bdy ! bdy e' il numero del body
	real time
 
        real deltat ! frazione di time step body si muove in ie
 
	integer ieltv(3,1)
        common /ieltv/ieltv !punta ad elemento opposto al nodo 

	real xgv(1),ygv(1) ! coordinate nodi
        integer nen3v(3,1) ! puntatore sui nodi di ogni elemento
        common /xgv/xgv, /ygv/ygv
        common /nen3v/nen3v

        real ulnv(nlvdim,1),vlnv(nlvdim,1)
        common /ulnv/ulnv, /vlnv/vlnv

	real v_ent ! valore mediato tra velocita' int e out
	real v_int ! modulo della velocita' di entrata
	real v_out ! modulo della velocita' di uscita
	integer l_int !lato di entrata numerazione elemento (var. output)
	integer l_out ! lato di uscita numerazione elemento
	real xbdy,ybdy ! coordinate del body
	real nxbdy,nybdy !nuove coordinate del body
	integer newie ! nuovo elemento che contiene il body

        integer idum ! seme f(x) random
                        
c variabili di servizio
	
	integer pb(3),ip,i,exi(2),ipb(2),exio(2)
	real nwdist,distance ! funzione per il calcolo del nwdist
        real dxx,dyy,ddxx,ddyy,ddd
        real itrx,itry
        real drd(3)
        real nni
        real xit,yit ! coordinate intercetta retta traj e retta enter
	real fsum,fpb(2),ffpb(2)
	real tcomp ! funzione per il calcolo del deltat
  	real ib,aax,icx,icy ! parametri retta del lato da cui entra il body
	real ocy,ocx,ob ! parametri retta del lato da cui esce il body
        real fay,fax ! parametri fascio rette traiettoria
	real db ! parametri retta della traiettoria sguita dal body 
	real dstbdy ! distanza tra lato entrata e uscita lungo la traiettoria
        real bdx ! distanza (frazione) del body da estremo piu vicno retta entrata
        real gf !precision along the side
        integer pbdx,gbdx ! numero interno dell'estremo piu vicino e piu lontano dal body
        integer ps,ng,nl,nnm,i1,i2
	data nnm/0/
	data pb/3*0/
        data idum/76342/
        save nnm
        save idum

c inizializzazione parametri

110 	continue
	nnm=nnm+1
        fsum=0

c==========================================================
c INDIVIDUAZIONE CASI
c==========================================================
        ps=0
        ng=0
        nl=0
        do i=1,3
         nni=vel_ie(i,ie)
         if(nni.lt.0)then
          ng=ng+1
          l_out=i
         elseif(nni.gt.0)then
          ps=ps+1
          l_int=i
         elseif(nni.eq.0)then
          nl=nl+1
         endif          
        end do

        v_ent=sqrt((ulnv(1,ie)**2)+(vlnv(1,ie)**2))
        
c_________________________________________________________
c CASO [A]: 1 flusso >0 (entrante), 2 flussi <0 (uscenti).
c Il lato di entrata � quello il cui flusso >0.
c - individuare fascio di rette traiettoria
c - individuare traiettoria del body
c - individuare lato di uscita
c - individuare nuove coordinate del body

        if((ng.eq.2).and.(ps.eq.1))then
         v_int=vel_ie(l_int,ie)
         
         if(v_int.le.0)then
          PRINT*,'STOP!! ERRORE VELOCITA DI ENTRATA NEGATIVA'
	  PRINT*,'CASO [A] TRACK_ORIG'
	  !stop
          ie = -ie
	  time=0
	  return
         endif                               
         
         
c calcolo retta del lato da cui � entrato il body
         
         i1=mod(l_int,3)+1
         i2=mod(i1,3)+1
	 exi(1)=i1
	 exi(2)=i2	 	
                                          
         
         call retta(exi,icy,icx,ib,ie)

         
c determinazione dei lati e dei flussi di uscita probabili pb(1,2),fpb(1,2) 
c pb(mi dice il numero interno del lato)


         ip=1
         do i=1,3
          if(l_int.ne.i)then
	   pb(i)=ip
           ipb(ip)=i
           fpb(ip)= vel_ie(i,ie)
           ip=ip+1
          endif
         enddo
         
         
c si calcola la frazione sul lato di entrata di competenza
c del lato pb(1) e del lato pb(2) per la possibile uscita del body         

        fsum=abs(fpb(1))+abs(fpb(2))
        do i=1,2
         ffpb(i)=abs(fpb(i)/(fsum))
        end do
       
       
c calcolo del fascio di rette traiettorie del body, mi serve sapere:
c (l_int) nodo opposto al lato di entrata coordinate del punto di 
c separazione (sep) tra le 2 frazioni sul lato di entrata (ad estremo 1 sommo ffpd(2))

        call dirz(l_int,ffpb,exi,icy,icx,ib,fay,fax,ie,pb,ipb)

        
c determinazione retta traiettoria
         
        call traj(fay,fax,xbdy,ybdy,db)


c deteminazione coordinate punto intercetto retta traiettoria-
c retta di entrata

        call interc(fay,fax,db,icy,icx,ib,xit,yit)


c calcola la distanza (bdx) del punto (xit,yit) dall'estremo piu vicino
c sulla retta di entrata in termini di frazione della lunghezza totale e
c dice qual'e' il numero interno dell'estremo piu vicino (pbdx) e 
c quello piu lontano (gbdx)

        call distp(xit,yit,exi,bdx,pbdx,gbdx,ie)
       
        if(ie.le.0)then
         call find_elem_from_old(-ie,xbdy,ybdy,ie)
         time=0
         return
        end if 
c individuo il lato da cui esce il body (l_out) mediante confronto
c tra la posizione xit,yit  e le distanze sulla retta di entrata
c individuate dalle frazioni (ffpd) 
        
        if(bdx.le.ffpb(pb(pbdx)))then
         l_out=pbdx
        elseif(bdx.gt.ffpb(pb(pbdx)))then
         l_out=gbdx
        end if

        
c determinazione velocit� di entrata v_ent del body v_ent
                
        v_out=vel_ie(l_out,ie)

        if(v_out.gt.0)then
         PRINT*,'STOP!! ERRORE VELOCIT� DI USCITA POSITIVA'
	 PRINT*,'CASO [A] TRACK_ORIG'	 
	 !stop
	 time=0
         ie = -ie
	 return
        elseif(v_out.eq.0)then
         PRINT*,'STOP!! ERRORE VELOCIT� DI USCITA NULLA'
	 PRINT*,'CASO [A] TRACK_ORIG'
	 !stop
	 time=0
         ie = -ie
	 return
        endif

        

c estremi del lato da cui esce il body

        i1=mod(l_out,3)+1
        i2=mod(i1,3)+1

        exio(1)=i1
        exio(2)=i2
        
        call retta(exio,ocy,ocx,ob,ie)

        
c calcolo punto intercetto tra retta traiettoria - retta uscita 

        call interc(fay,fax,db,ocy,ocx,ob,itrx,itry)

c___________________________________________________________
c CASO [B]: 2 flussi >0 (entranti) o 1 flusso >0
c 1 flusso =0, 1 flusso <0 (uscente).
c Il lato di uscita � quello il cui flusso <0.
c - individuare fascio di rette traiettoria
c - individuare traiettoria del body
c - individuare lato di entrata
c - individuare nuove coordinate del body

       elseif(((ng.eq.1).and.(ps.eq.2)).or.
     +((ng.eq.1).and.(nl.eq.1).and.(ps.eq.1)))then
        v_out=vel_ie(l_out,ie)

        
        if(v_out.ge.0)then
         PRINT*,'STOP!! ERRORE VELOCITA DI USCITA POSITIVO'
         PRINT*,'CASO [B] TRACK_ORIG'
	 !stop
	 time=0
         ie = -ie
	 return
        endif

c calcolo retta del lato da cui uscir� il body

	i1=mod(l_out,3)+1
        i2=mod(i1,3)+1

                
        exio(1)=i1
        exio(2)=i2
        
        call retta(exio,ocy,ocx,ob,ie)

        
c determinazione dei lati e dei flussi di entrata probabili pb(1,2),
c fpb(1,2), pb(mi dice il numero interno del lato)
        
        ip=1
        do i=1,3
         if(l_out.ne.i)then
          pb(i)=ip
          ipb(ip)=i
          fpb(ip)=vel_ie(i,ie)
          ip=ip+1
         endif
        end do
        
        
c si calcola la frazione sul lato di uscita, di competenza
c del lato pb(1) e del lato pb(2) per la possibile entrata del body
                             
        fsum=abs(fpb(1))+abs(fpb(2))
        do i=1,2
         ffpb(i)=abs(fpb(i)/fsum)
        end do
     
        
c calcolo del fascio di rette traiettorie del body, mi serve sapere:
c (l_out) nodo opposto al lato di uscita coordinate del punto di
c separazione (sep) tra le 2 frazioni sul lato di entrata (ad estremo 1 sommo ffpd(2))

        call dirz(l_out,ffpb,exio,ocy,ocx,ob,fay,fax,ie,pb,ipb)
        
c determinazione retta traiettoria

        call traj(fay,fax,xbdy,ybdy,db)

c deteminazione coordinate punto intercetto retta traiettoria-
c retta di uscita

        call interc(fay,fax,db,ocy,ocx,ob,itrx,itry)
        
c calcola la distanza (bdx) del punto (itrx,itry) dall'estremo piu vicino
c sulla retta di uscita in termini di frazione della lunghezza totale e
c dice qual'e' il numero interno dell'estremo piu vicino (pbdx) e
c quello piu lontano (gbdx)
        
        call distp(itrx,itry,exio,bdx,pbdx,gbdx,ie)

        if(ie.le.0)then
         call find_elem_from_old(-ie,xbdy,ybdy,ie)
         time=0
         return
        end if
                                                           
c individuo il lato da cui entra il body (l_int) mediante confronto
c tra la posizione itrx, itry  e le distanze sulla retta di uscita
c individuate dalle frazioni (ffpd)

        if(bdx.le.ffpb(pb(pbdx)))then
         l_int=pbdx
        elseif(bdx.gt.ffpb(pb(pbdx)))then
         l_int=gbdx
        end if


c determinazione velocit� di entrata v_ent del body v_ent
         
        v_int=vel_ie(l_int,ie)
                
        if(v_int.lt.0)then
         PRINT*,'STOP!! ERRORE VELOCIT� DI ENTRATA NEGATIVA'
         PRINT*,'CASO [B] TRACK_ORIG'
	 !stop
	 time=0
         ie = -ie
	 return
        elseif(v_int.eq.0)then
         PRINT*,'STOP!! ERRORE VELOCIT� DI ENTRATA NULLA'
	 PRINT*,'CASO [B] TRACK_ORIG'
	 !stop
	 time=0
         ie = -ie
	 return
        endif
c____________________________________________________________________________
c CASO [C]:tutti i flussi sono entranti >0 O NULLI =0
c Errore il campo di moto � convergente, problemi con l'eliminazione delle componenti
c DIVERGENTI FERMARE IL CILO

       elseif((ps.eq.3).or.(nl.eq.3))then
c       print*,'STOP!! FLUSSI TUTTI POSITIVI O NULLI IN ',ie
c       PRINT*,'CASO [C] TRACK_ORIG'
         time=0
c        ie = -ie
        return
       elseif((ng.eq.3))then
c       print*,'STOP!! FLUSSI TUTTI NEGATIVI IN ',ie
c       PRINT*,'CASO [D] TRACK_ORIG'
         time=0
c        ie = -ie
        return
       elseif((nl.eq.2).and.(ps.eq.1))then
c       print*,'STOP!! UN SOLO FLUSSO POSITIVO E 2 NULLI IN ',ie
c       PRINT*,'CASO [E] TRACK_ORIG'
c       PRINT*,'[elemento di partenza con due lati C.B]'
         time=0
c        ie = -ie
        return
       elseif((nl.eq.2).and.(ng.eq.1))then
c       print*,'STOP!! UN SOLO FLUSSO NEGATIVO E 2 NULLI IN ',ie
c       PRINT*,'CASO [F] TRACK_ORIG'
c       PRINT*,'[elemento di partenza con due lati C.B]'
         time=0
c        ie = -ie
        return
       elseif((nl.eq.1).and.(ps.eq.2))then
c       print*,'STOP!! 2 FLUSSI POSITIVI E 1 NULLO IN ',ie
c       PRINT*,'CASO [G] TRACK_ORIG'
         time=0
c        ie = -ie
        return
       elseif((nl.eq.1).and.(ng.eq.2))then
c       print*,'STOP!! 2 FLUSSI NEGATIVI E 1 NULLO IN ',ie
c       PRINT*,'CASO [H] TRACK_ORIG'
         time=0
c        ie = -ie
        return
       end if         
c_________________________________________________________________________


c calcolo distanza massima percorribile all'interno di elemento
       
       dxx=itrx-xbdy
       dyy=itry-ybdy
       ddxx=(dxx)**2
       ddyy=(dyy)**2
       ddd=ddxx+ddyy
       dstbdy=sqrt(ddd)
       nwdist=time*v_ent 
        
c=====================================================================
c CALCOLO NUOVE COORDINATE DEL BODY (nxbdy,nybdy)
c=====================================================================

c calcolo della distanza percorsa in 1 tstep (nwdist)
c se il body arriva esattamente sul lato a ttime=0 (1 iflogico)?
c se si esce dall'elemento (2 iflogico) allora ho coordinate del body sul lato di uscita
c e in piu un deltat > 0 da utilizzare durnate questo stesso tstep.
c se si rimane nell'elemento (3 iflogico) allora calcolo le nuove coordinate 
c (pnt_inside) e passo al tstep successivo

c Se il body finisce a 0.1 m dalla linea decido si farlo passare al
c elemento ad una distanza di 0.1 m. Questa � la precisione da me 
c imposta
        gf=dvert(1,ie)/100000 
        if(abs(nwdist-dstbdy).lt.gf)then      ! body dista meno di 0.1m dal lato
          PRINT*,'WARNING!! BODY SU LATO ELEMENTO'
          PRINT*,'track orig',bdy
          nwdist=dstbdy+2*gf
        endif
        if(nwdist.gt.dstbdy)then            !2 body su nuovo el
	 nxbdy=itrx
	 nybdy=itry
	 distance=nwdist-dstbdy
	 deltat=distance/v_ent
         newie=ieltv(l_out,ie)
         if(newie.eq.-1)then			!2.a body e uscito dal dominio
	  if( .not. bback ) then
            PRINT*,'STOP!! BODY ',bdy,' USCITO '
            PRINT*,'ELEMENTO USCITA ',ie
	  end if
          ie=-ie ! flag per skip bdy dal calcolo
	  time=0.
          xbdy=nxbdy ! ultime coordinate del body uscito
          ybdy=nybdy ! ultime coordinate del body uscito
          return
         endif
         do i=1,3
          if(ie.eq.ieltv(i,newie))then         
           lt_body(bdy)=i ! individuazione lato di entrata del body
          endif
         end do
         xbdy=nxbdy ! nuove coordinate del body
	 ybdy=nybdy ! nuove coordinate del body
         ie=newie ! individuazione nuovo elemento
         time=deltat
  	elseif(nwdist.lt.dstbdy)then            !3 body su vecchio el
	 call pnt_inside(nwdist,dstbdy,fay,fax,xbdy,ybdy
     +				,itrx,itry,nxbdy,nybdy)
	 deltat=0.
	 time=deltat
	 ie=ie ! individuazione vecchio elemento
         xbdy=nxbdy ! nuove coordinate del body
         ybdy=nybdy ! nuove coordinate del body
        endif
        end

c**********************************************************************

	subroutine track_line(time,bdy,ie,xbdy,ybdy)

c in questa subroutine si calcola il percorso del body partendo
c da un punto sul lato dell'elemento ie
c questa subroutine e' chiamata ogni qualvolta un body
c nello stesso timestep arriva in un nuovo elemento

        implicit none
        
	include 'param.h' 
	include 'lagrange.h' 

	integer ie,bdy ! bdy e' il numero del body
	real time
 
        real deltat ! frazione di time step body si muove in ie
 
	integer ieltv(3,1)
        common /ieltv/ieltv !punta ad elemento opposto al nodo 

	real xgv(1),ygv(1) ! coordinate nodi
        integer nen3v(3,1) ! puntatore sui nodi di ogni elemento
        common /xgv/xgv, /ygv/ygv
        common /nen3v/nen3v

        real ulnv(nlvdim,1),vlnv(nlvdim,1)
        common /ulnv/ulnv, /vlnv/vlnv
                
	real v_ent ! valore mediato tra velocita' int e out
	real v_int ! modulo della velocita' di entrata
	real v_out ! modulo della velocita' di uscita
	integer l_int !lato di entrata numerazione elemento (var. input)
	integer l_out ! lato di uscita numerazione elemento
	real xbdy,ybdy ! coordinate del body
	real nxbdy,nybdy !nuove coordinate del body
	integer newie ! nuovo elemento che contiene il body

        integer idum ! seme f(x) random
                        

c variabili di servizio
	
	integer pb(3),ip,i,exi(2),ipb(2),exio(2)
	real nwdist,distance ! funzione per il calcolo del nwdist
        real dxx,dyy,ddxx,ddyy,ddd
        real itrx,itry
        real drd(3)
        real fsum,fpb(2),ffpb(2)
	real tcomp ! funzione per il calcolo del deltat
  	real icy,icx,ib ! parametri retta del lato da cui entra il body
	real ocy,ocx,ob ! parametri retta del lato da cui esce il body
        real fay,fax ! parametri fascio rette traiettoria
        real db ! parametri retta della traiettoria sguita dal body 
	real dstbdy ! distanza tra lato entrata e uscita lungo la traiettoria
        real bdx ! distanza (frazione) del body da estremo piu vicno retta entrata
	integer pbdx,gbdx ! numero interno dell'estremo piu vicino e piu lontano dal body
        integer ps,ng,nl,i1,i2
        real gf !precision along the side
c inizializzazione parametri
        data idum/947437/
        save idum
        
110 	continue

        fsum=0
        
c==========================================================
c CALCOLO TRAIETTORIE
c==========================================================

c calcolo retta del lato da cui entra il body
c exi(1,2) sono i puntatori del nen3v per gli estremi della retta 
c su cui si trova il body

        
        l_int=lt_body(bdy)


        i1=mod(l_int,3)+1
        i2=mod(i1,3)+1
        exi(1)=i1
        exi(2)=i2
                                          
	call retta(exi,icy,icx,ib,ie)

c velocit� di entrata v_int 

	v_int=vel_ie(l_int,ie)

        if(v_int.le.0)then
         PRINT*,'STOP!! ERRORE FLUSSO DI ENTRATA NEGATIVO'
         print*,'track line'
	 !stop
	 time=0
         ie = -ie
	 return
        endif

c determinazione dei lati e dei flussi di uscita probabili pb(1,2),vpb(1,2)
c pb(mi dice il numero interno del lato)
        
        ip=1
        do i=1,3
         if(l_int.ne.i)then
 	  pb(i)=ip
	  ipb(ip)=i
          fpb(ip)=vel_ie(i,ie) 
          ip=ip+1  
         endif 
	end do

        v_ent=sqrt((ulnv(1,ie)**2)+(vlnv(1,ie)**2))
         
         
c______________________________________________________________________
c CASO [A] flussi ai lati opposti al l_int entrambi <0
c si calcola la frazione sul lato di entrata di competenza 
c del lato pb(1) e del lato pb(2) per la possibile uscita del body.

      	if((fpb(1).lt.0).and.(fpb(2).lt.0))then

                
        fsum=abs(fpb(1))+abs(fpb(2))
         do i=1,2
	  ffpb(i)=abs(fpb(i)/fsum)
	 end do

         
c calcola la distanza del body (bdx) dall'estremo piu vicino
c in termini di frazione sulla retta totale  e
c dice qual'e' il numero interno dell'estremo piu vicino (pbdx)
c e quello piu lontano (gbdx)

        call distp(xbdy,ybdy,exi,bdx,pbdx,gbdx,ie)
        if(ie.le.0)then
         call find_elem_from_old(-ie,xbdy,ybdy,ie)
         time=0
         return
        end if
                                                   
        
c individuo il lato da cui esce il body mediante confronto
c tra posizione del body rispetto alle frazioni  

	 if(bdx.le.ffpb(pb(pbdx)))then	  
	  l_out=pbdx
	 elseif(bdx.gt.ffpb(pb(pbdx)))then
	  l_out=gbdx
	 end if

         
c estremi del lato da cui esce il body

        i1=mod(l_out,3)+1
        i2=mod(i1,3)+1

        exio(1)=i1
        exio(2)=i2

                                          
c equazione della retta del lato di uscita

        call retta(exio,ocy,ocx,ob,ie)


c calcolo del fascio di traiettorie del body, mi serve sapere:
c (l_int) nodo opposto al lato di entrata
c coordinate del punto di separazione (sep) tra le 2 frazioni sul lato
c di entrata (ad estremo 1 sommo ffpd(2))
         
        call dirz(l_int,ffpb,exi,icy,icx,ib,fay,fax,ie,pb,ipb)
        
c determinazione retta traiettoria

        call traj(fay,fax,xbdy,ybdy,db)

c deteminazione coordinate punto intercetto retta traiettoria-
c retta di uscita
         
        call interc(fay,fax,db,ocy,ocx,ob,itrx,itry)

c_____________________________________________________________________________
c CASO [B] flussi ai lati opposti al l_int uno >0 e altro <0
c in questo caso il lato di uscita e' sicuramente quello dove ho la velocita'
c negativa, non serve quindi il calcolo delle rette di partenza e di arrivo
c gli estremi della retta del lato di uscita sono i due estremi non opposti
c alla retta che ha valore di vel negativo

        else if(((fpb(1).lt.0).and.(fpb(2).ge.0)).or.
     + ((fpb(1).ge.0).and.(fpb(2).lt.0)))then
        

         do i=1,2
          if(fpb(i).lt.0)then
           l_out=ipb(i)
          endif
         enddo

c estremi del lato da cui esce il body
	
        i1=mod(l_out,3)+1
        i2=mod(i1,3)+1

        exio(1)=i1
        exio(2)=i2
                                          
                                          
c equazione della retta del lato di uscita

        call retta(exio,ocy,ocx,ob,ie)
        
        
c individuo i 2 flussi di entrata del body nell'elemento per il calcolo
c delle frazioni di competenza sul lato di uscita.
 
         ip=1
         do i=1,3
          if(l_out.ne.i)then
           pb(i)=ip
           ipb(ip)=i
           fpb(ip)=vel_ie(i,ie)
           ip=ip+1
          endif
         end do

         
c calcolo delle frazioni sul lato di uscita per determinare 
c la retta traiettoria del percorso del body
	
         fsum=abs(fpb(1))+abs(fpb(2))
         do i=1,2
          ffpb(i)=abs(fpb(i)/fsum)
         end do


c calcolo del fascio di traiettorie del body, mi serve sapere:
c l_out nodo opposto al lato di uscita
c coordinate del punto di separazione tra le 2 frazioni sul lato 
c di uscita (come prima ad estremo 1 sommo ffpb(2))

         call dirz(l_out,ffpb,exio,ocy,ocx,ob,fay,fax,ie,pb,ipb)
         
c calcolo retta traiettoria 

         call traj(fay,fax,xbdy,ybdy,db)
         
c deteminazione coordinate punto intercetto retta traiettoria-
c retta di uscita

         call interc(fay,fax,db,ocy,ocx,ob,itrx,itry)

c calcolo velocita' in uscita v_out per il calcolo del modulo della
c velocita di percorso v_ent pari alla media v_int, v_out

         v_out=vel_ie(l_out,ie)
        
         if(v_out.gt.0)then
          PRINT*,'STOP!! ERRORE VELOCIT� DI USCITA POSITIVA'
          PRINT*,'track line'
	  !stop
	  time=0
          ie = -ie
	  return
         elseif(v_out.eq.0)then
          PRINT*,'STOP!! ERRORE VELOCIT� DI USCITA NULLA'
          PRINT*,'track line'
	  !stop
	  time=0
          ie = -ie
	  return
         endif

        elseif((fpb(1).ge.0).and.(fpb(2).ge.0))then
c________________________________________________________________________________
c CASO [C] tutti i flussi sono entranti o nulli. Errore il campo di moto �
c convergente, problemi con l'eliminazione delle componenti DIVERGENTI
c FERMARE IL CILO

c        print*,'STOP!! FLUSSI TUTTI POSITIVI IN ',ie
c        PRINT*,'CASO [C] TRACK_LINE'
         time=0
c        ie = -ie
         return
        elseif((fpb(1).lt.0).and.(fpb(2).lt.0))then
c        print*,'STOP!! FLUSSI TUTTI NEGATIVI IN ',ie
c        PRINT*,'CASO [D] TRACK_LINE'
         time=0
c        ie = -ie
         return
        elseif((nl.eq.2).and.(ps.eq.1))then
c        print*,'STOP!! UN SOLO FLUSSO POSITIVO E 2 NULLI IN ',ie
c        PRINT*,'CASO [E] TRACK_LINE'
c        PRINT*,'[elemento con due lati C.B]'
         time=0
c        ie = -ie
         return
        elseif((nl.eq.2).and.(ng.eq.1))then
c        print*,'STOP!! UN SOLO FLUSSO NEGATIVO E 2 NULLI IN ',ie
c        PRINT*,'CASO [F] TRACK_LINE'
c        PRINT*,'[elemento con due lati C.B]'
         time=0
c        ie = -ie
         return
        elseif((nl.eq.1).and.(ps.eq.2))then
c        print*,'STOP!! 2 FLUSSI POSITIVI E 1 NULLO IN ',ie
c        PRINT*,'CASO [G] TRACK_LINE'
         time=0
c        ie = -ie
         return
        elseif((nl.eq.1).and.(ng.eq.2))then
c        print*,'STOP!! 2 FLUSSI NEGATIVI E 1 NULLO IN ',ie
c        PRINT*,'CASO [H] TRACK_LINE'
         time=0
c        ie = -ie
         return
        endif
c_______________________________________________________________________

c calcolo distanza massima percorribile all'interno di elemento

        dxx=itrx-xbdy
        dyy=itry-ybdy
        ddxx=(dxx)**2
        ddyy=(dyy)**2
        ddd=ddxx+ddyy
        dstbdy=sqrt(ddd)

	nwdist=time*v_ent


c=====================================================================
c CALCOLO NUOVE COORDINATE DEL BODY (nxbdy,nybdy)
c=====================================================================

c il body entra sempre da un lato con una velocita' positiva
c calcolo della distanza percorsa in 1 tstep (nwdist)
c se il body arriva esattamente sul lato a ttime=0 (1 iflogico)?
c se si esce dall'elemento (2 iflogico) allora ho coordinate del body sul lato di uscita
c e in piu un deltat > 0 da utilizzare durnate questo stesso tstep.
c se si rimane nell'elemento (3 iflogico) allora calcolo le nuove coordinate 
c (pnt_inside) e passo al tstep successivo

c Se il body finisce a 0.1 m dalla linea decido si farlo passare al
c elemento ad una distanza di 0.1 m. Questa � la precisione da me
c imposta

        gf=dvert(1,ie)/100000    
        if(abs(nwdist-dstbdy).lt.gf)then     ! body dista meno di 0.1m dal lato    
          PRINT*,'STOP!! BODY SU LATO ELEMENTO'
          PRINT*,'track line',bdy
          nwdist=dstbdy+2*gf
        endif

        if(nwdist.gt.dstbdy)then            !2 body su nuovo el
	 nxbdy=itrx
	 nybdy=itry
	 distance=nwdist-dstbdy
         deltat=distance/v_ent
         newie=ieltv(l_out,ie) !
         if(newie.eq.-1)then			!2.a body uscito dal dominio
	  if( .not. bback ) then
            PRINT*,'STOP!! BODY ',bdy,' USCITO '
            PRINT*,'ELEMENTO USCITA ',ie
	  end if
          ie=-ie! flag per skip bdy dal calcolo
          xbdy=nxbdy ! ultime coordinate del body uscito
          ybdy=nybdy ! ultime coordinate del body uscito
	  time=0.
	  return     
         endif
         do i=1,3
          if(ie.eq.ieltv(i,newie))then          !3 body su vecchio el
           lt_body(bdy)=i ! individuazione lato di entrata del body
          endif
         end do
         xbdy=nxbdy ! nuove coordinate del body
	 ybdy=nybdy ! nuove coordinate del body
        ie=newie ! individuazione nuovo elemento
        time=deltat
       elseif(nwdist.lt.dstbdy)then
	 call pnt_inside(nwdist,dstbdy,fay,fax,xbdy,ybdy
     +				,itrx,itry,nxbdy,nybdy)
	 deltat=0.
	 time=deltat
	 ie=ie ! individuazione vecchio elemento
         xbdy=nxbdy ! nuove coordinate del body
         ybdy=nybdy ! nuove coordinate del body
       endif
       end

c**********************************************************************

