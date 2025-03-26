seed_summary = list(
  # SIMULATION_1
  TVpsi1_pos = seed_nul.norm5,
  TVpsi1_neg = seed_pos.psi1,
  
  SIM1_pos = seed_neg.psi1,
  SIM1_neg = seed_DR.pos2.norm,
  SIM1_nul = seed_DR.neg2.norm,
  
  # SIMULATION_2
  SIM2_neg1 = seed_DR.nul.norm,
  SIM2_neg2 = seed_neg.norm1,
  SIM2_neg3 = seed_neg.norm2,
  SIM2_neg4 = seed_neg.norm3,
  SIM2_neg5 = seed_neg.norm4,
  
  SIM2_nul1 = seed_neg.norm5,
  SIM2_nul2 = seed_nul.norm1,
  SIM2_nul3 = seed_nul.norm2,
  SIM2_nul4 = seed_nul.norm3,
  SIM2_nul5 = seed_nul.norm4,
  
  # SIMULATION_1_BOUNDED
  SIMS1_pos = SIM.RES_DR$seed_ALT.BOTHmis,
  SIMS1_neg = SIM.RES1_DR$seed_pos,
  SIMS1_nul = SIM.RES1_DR$seed_neg,
  
  # SIMULATION_2_BOUNDED
  SIMS2_pos1 = SIM.RES1_DR$seed_nul,
  SIMS2_pos2 = seed_pos1,
  SIMS2_pos3 = seed_pos2,
  SIMS2_pos4 = seed_pos3,
  
  SIMS2_neg1 = seed_pos4,
  SIMS2_neg2 = seed_neg1,
  SIMS2_neg3 = seed_neg2,
  SIMS2_neg4 = seed_neg3,
  SIMS2_neg5 = seed_neg4,
  
  SIMS2_nul1 = seed_neg5,
  SIMS2_nul2 = seed_nul1,
  SIMS2_nul3 = seed_nul2,
  SIMS2_nul4 = seed_nul3,
  SIMS2_nul5 = seed_nul4,
  
  # SIMULATION_S7
  SIMS7_discrete = seed_DR.nul.norm_half,
  SIMS7_continuous = seed_Wang
)

save(list = "seed_summary",
     file = paste0("./SubmitPrograms/seed_summary.RData"))
