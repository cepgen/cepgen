      subroutine grv95lo(xpart,q2part,uv,dv,us,ds,ss,wg)

      amulo=0.23
      argulo=0.232*0.232
      rat1 = q2part/argulo
      rat2 = amulo/argulo
      s=alog(alog(rat1)/alog(rat2))
      s2=s*s
      s3=s2*s
      s12=sqrt(s)
      x1=1.-xpart
      x12=sqrt(xpart)
      x32=x12*xpart
      alx=alog(1./xpart)
c
c-----uvalence
c
      auv=0.59-0.024*s
      buv=0.131+0.063*s
      anuv=2.284+0.802*s+0.055*s2
      aduv=-0.449-0.138*s-0.076*s2
      bduv=0.213+2.669*s-0.728*s2
      cduv=8.854-9.135*s+1.979*s2
      dduv=2.997+0.753*s-0.076*s2
c-------
      ucoff=anuv*xpart**auv
      uwna=1.+aduv*xpart**buv+bduv*xpart+cduv*x32
      xuv1=x1**dduv
      uv=ucoff*uwna*xuv1
c
c-----dvalence
c
      adv=0.376
      bdv=0.486+0.062*s
      andv=0.371+0.083*s+0.039*s2
      addv=-0.509+3.31*s-1.248*s2
      bddv=12.41-10.52*s+2.267*s2
      cddv=6.373-6.208*s+1.418*s2
      dddv=3.691+0.799*s-0.071*s2
c-----
      dcoff=andv*xpart**adv
      dwna=1.+addv*xpart**bdv+bddv*xpart+cddv*x32
      xdv1=x1**dddv
      dv=dcoff*dwna*xdv1
c
c-----asymmetric sea
c
      adel=0.409-0.005*s
      bdel=0.799+0.071*s
      andel=0.082+0.014*s+0.008*s2
      addel=-38.07+36.13*s-0.656*s2
      bddel=90.31-74.15*s+7.645*s2
      cddel=0.
      dddel=7.486+1.217*s-0.159*s2
c----------------------------------
      asymc=andel*xpart**adel
      awndel=1.+addel*xpart**bdel+bddel*xpart
      awndel=awndel+cddel*x32
      del=asymc*awndel*x1**dddel
c
c-----symmetric sea
c
      alphs=1.451
      bets=0.271
      as=0.41-0.232*s
      bs=0.534-0.457*s
      ads=0.89-0.14*s
      bds=-0.981
      cds=0.32+0.683*s
      dds=4.752+1.164*s+0.286*s2
      eds=4.119+1.713*s
      edprs=0.682+2.978*s
c -----
      awn1=xpart**as
      awn2=ads+bds*xpart+cds*xpart*xpart
      aw1=awn1*awn2*alx**bs
      awn3=exp(-eds+sqrt(edprs*s**bets*alx))
      ws=(aw1+s**alphs*awn3)*x1**dds
c------
c-----strange sea
c------
      alphss=0.914
      betss=0.577
      ass=1.798-0.596*s
      adss=-5.548+3.669*s12-0.616*s
      bdss=18.92-16.73*s12+5.168*s
      ddss=6.379-0.35*s+0.142*s2
      edss=3.981+1.638*s
      edprss=6.402
c---------------------------------
      ssc=s**alphss/alx**ass
      awn=1.+adss*x12+bdss*xpart
      sdlo=exp(-edss+sqrt(edprss*s**betss*alx))
      ss=ssc*awn*x1**ddss*sdlo
c
c-----note that strange sea=ss=s+sbar
c
c-----nonstrange sea
c
      us=(ws-del)/2.
      ds=(ws+del)/2.
c
c-----gluons
c
      alphg=0.524
      betg=1.088
      ag=1.742-0.930*s
      bg=-0.399*s2
      adg=7.486-2.185*s
      bdg=16.69-22.74*s+5.779*s2
      cdg=-25.59+29.71*s-7.296*s2
      ddg=2.792+2.215*s+0.422*s2-0.104*s3
      edg=0.807+2.005*s
      edprg=3.841+0.316*s
c -----
      awn1=xpart**ag
      awn2=adg+bdg*xpart+cdg*xpart*xpart
      aw1=awn1*awn2*alx**bg
      awn3=exp(-edg+sqrt(edprg*s**betg*alx))
      wg=(aw1+s**alphg*awn3)*x1**ddg

      return
      end
