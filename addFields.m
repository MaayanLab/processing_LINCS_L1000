function [chdir] = addFields(chdir,proto)
    chdir.batch = proto.batch;
    chdir.pert_id = proto.SM_LINCS_ID;
    chdir.pert_dose = proto.SM_Dose;
    chdir.sig_id = strjoin({proto.batch,proto.SM_LINCS_ID,proto.SM_Dose},':');
    
    chdir.pert_desc = proto.SM_Name;
    chdir.pert_type = proto.SM_Pert_Type;
    chdir.pert_time = proto.SM_Time;
    chdir.pert_time_unit = proto.SM_Time_Unit;
    chdir.pert_dose_unit = proto.SM_Dose_Unit;
    chdir.cell_id = proto.CL_Name;

   
end