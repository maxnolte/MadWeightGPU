ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DOUBLE COMPLEX FUNCTION COND(CONDITION,TRUECASE,FALSECASE)
      IMPLICIT NONE
      DOUBLE COMPLEX CONDITION,TRUECASE,FALSECASE
      IF(CONDITION.EQ.(0.0D0,0.0D0)) THEN
        COND=TRUECASE
      ELSE
        COND=FALSECASE
      ENDIF
      END

      DOUBLE COMPLEX FUNCTION REGLOG(ARG)
      IMPLICIT NONE
      DOUBLE COMPLEX ARG
      IF(ARG.EQ.(0.0D0,0.0D0)) THEN
        REGLOG=(0.0D0,0.0D0)
      ELSE
        REGLOG=LOG(ARG)
      ENDIF
      END
