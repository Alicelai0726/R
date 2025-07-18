# Survival Modeling in Oncology Using R
A comprehensive survival analysis and economic evaluation for a hypothetical cancer treatment, Ultradrug, using R. 
The analysis follows the structure of a health technology assessment, incorporating survival modeling, treatment switching adjustments, and cost-effectiveness modeling. # Non-parametric methods such as Kaplan–Meier curves and RMST were used to compare treatment arms. 
Proportional hazard assumptions were assessed using diagnostic plots including log–log plots, Schoenfeld residuals, and smoothed hazard functions. Parametric models—both combined and independently fitted—were evaluated using the flexsurv package, with log-normal selected as the base-case model due to its flexible hazard representation and strong visual fit. Additional packages such as survival, survminer, survRM2, bshazard, ggplot2, and muhaz supported data processing, visualization, and model diagnostics. 
Treatment switching was adjusted using a Two-Stage Estimation approach based on an AFT model, and a partitioned survival model was built to estimate QALYs and ICERs.

# Cost-Utility Analysis of Enhanced Physiotherapy for Motor Neuron Disease Using EQ-5D and SF-6D: A Missing Data and Perspective-Based Evaluation
This project conducts a comprehensive cost-utility analysis (CUA) of enhanced physiotherapy for motor neuron disease (MND) using R. 
The analysis began with data preparation and handling of substantial missingness in key variables (EQ-5D: 38%, SF-6D: 15%, costs: 30%). 
Multiple imputation by chained equations (MICE) was applied using predictive mean matching (PMM), with further comparisons using random forest and hot-deck methods to assess imputation sensitivity. 
The study evaluated QALYs based on EQ-5D and SF-6D at multiple time points (0, 3, 12, and 24 months) using the trapezium rule, with costs estimated under both NHS and hypothetical societal perspectives. Cost-effectiveness was assessed through ICERs, cost-effectiveness planes, and CEACs.

