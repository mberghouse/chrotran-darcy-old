# swap_comment_character.py
import sys
import shutil
import os
import fnmatch


def append_dict_to_list(list_,dict_):
    for key, value in dict_.items():
        if value:
            if value == True:
                list_.append('{}\n'.format(key))
            else:
                list_.append('{} {}\n'.format(key,value))
    
def write_list(f,list_):
    if len(list_) > 0:
        for line in list_:
            f.write('  '+line)
        f.write('\n') 
        
def write_block(f,process,list_):    
    if list_[0].strip() != process:
        list_.insert(0,'{}\n'.format(process))
    for i in range(1,len(list_)):
        list_[i] = '  ' + list_[i]
    list_.append('/\n')
    write_list(f,list_)

def refactor_file(filename,replace_file_flag):
    timestepper_flow_dict = {
        'MAX_PRESSURE_CHANGE' : None,
        'MAX_TEMPERATURE_CHANGE' : None,
        'MAX_CONCENTRATION_CHANGE' : None,
        'MAX_CFL' : None,
        'MAX_SATURATION_CHANGE' : None,
        'GAS_SAT_THRESH_FORCE_TS_CUT' : None,
        'MIN_LIQ_PRES_FORCE_TS_CUT' : None,
        'MAX_ALLOW_GAS_SAT_CHANGE_TS' : None,
        'MAX_ALLOW_LIQ_PRES_CHANGE_TS' : None,
        'GAS_SAT_CHANGE_TS_GOVERNOR' : None,
        'LIQ_PRES_CHANGE_TS_GOVERNOR' : None,
        'GAS_SAT_GOV_SWITCH_ABS_TO_REL' : None,
        'MINIMUM_TIMESTEP_SIZE' : None
        }

    timestepper_trans_dict = {
        'MAX_CFL' : None,
        'MAX_VOLUME_FRACTION_CHANGE' : None
        }

    newton_flow_dict = {
        'NUMERICAL_JACOBIAN' : None,
        'ANALYTICAL_JACOBIAN' : None,
        'MAX_NEWTON_ITERATIONS' : None,
        'RESIDUAL_INF_TOL' : None,
        'RESIDUAL_ABS_INF_TOL' : None,
        'LIQUID_RESIDUAL_ABS_INF_TOL' : None,
        'GAS_RESIDUAL_ABS_INF_TOL' : None,
        'ENERGY_RESIDUAL_ABS_INF_TOL' : None,
        'RESIDUAL_SCALED_INF_TOL' : None,
        'ITOL_SCALED_RESIDUAL' : None,
        'LIQUID_RESIDUAL_SCALED_INF_TOL' : None,
        'GAS_RESIDUAL_SCALED_INF_TOL' : None,
        'ENERGY_RESIDUAL_SCALED_INF_TOL' : None,
        'UPDATE_INF_TOL' : None,
        'ABS_UPDATE_INF_TOL' : None,
        'PRES_ABS_UPDATE_INF_TOL' : None,
        'TEMP_ABS_UPDATE_INF_TOL' : None,
        'SAT_ABS_UPDATE_INF_TOL' : None,
        'XMOL_ABS_UPDATE_INF_TOL' : None,
        'LIQUID_PRES_ABS_UPDATE_INF_TOL': None,
        'GAS_PRES_ABS_UPDATE_INF_TOL' : None,
        'AIR_PRES_REL_UPDATE_INF_TOL' : None,
        'PRESSURE_DAMPENING_FACTOR' : None,
        'MAX_ITERATION_BEFORE_DAMPING' : None,
        'DAMPING_FACTOR' : None,
        'ITOL_RELATIVE_UPDATE' : None,
        'PRESSURE_CHANGE_LIMIT' : None,
        'SATURATION_CHANGE_LIMIT' : None,
        'TEMPERATURE_CHANGE_LIMIT' : None,
        'USE_INFINITY_NORM_CONVERGENCE' : None,
        'LIQUID_RESIDUAL_INFINITY_TOL' : None,
        'GAS_RESIDUAL_INFINITY_TOL' : None,
        'MAX_ALLOW_REL_LIQ_PRES_CHANGE_NI' : None,
        'MAX_ALLOW_REL_GAS_SAT_CHANGE_NI' : None,
        'GAS_SAT_THRESH_FORCE_EXTRA_NI' : None
        }

    newton_trans_dict = {
        'NUMERICAL_JACOBIAN' : None,
        'ANALYTICAL_JACOBIAN' : None,
        'ITOL_RELATIVE_UPDATE' : None
        }

    f = open(filename,'r')
    
    linear_solver_flow = []
    newton_solver_flow = []
    timestepper_flow = []
    
    linear_solver_tran = []
    newton_solver_tran = []
    timestepper_tran = []
    
    tupl = ('TIMESTEPPER','NEWTON_SOLVER','LINEAR_SOLVER')
    
    
    
#    print('here1')
    store_mode = False
    while True:
        line = f.readline()
#        print(line)
#        print(len(line))
        if len(line) == 0:
            break
        w = line.split()
        card = ''
        if len(w) > 0:
          card = w[0].strip().upper()

        if card in timestepper_flow_dict:
            if len(w) > 1:
                timestepper_flow_dict[card] = w[1]
            else:
                timestepper_flow_dict[card] = True
        if card in timestepper_trans_dict:
            if len(w) > 1:
                timestepper_trans_dict[card] = w[1]
            else:
                timestepper_trans_dict[card] = True
        if card in newton_flow_dict:
            if len(w) > 1:
                newton_flow_dict[card] = w[1]
            else:
                newton_flow_dict[card] = True
        if card in newton_trans_dict:
            if len(w) > 1:
                newton_trans_dict[card] = w[1]
            else:
                newton_trans_dict[card] = True
        if card.startswith(tupl):
            store_mode = True
            list_ = []
            block_card = card
            
            if len(w) > 1:
                card2 = w[1].strip().upper()
            else:
                card2 = ""
            if card2 == 'FLOW':
                process_model = 1
                line = '{} \n'.format(card)
            elif card2 == 'TRANSPORT':
                process_model = 2
                line = '{} \n'.format(card)
            else:
                process_model = 0
        
        elif store_mode and \
             (card.startswith('END') or card.startswith('/')):
            store_mode = False
            if block_card.startswith('TIMESTEPPER'):
                if process_model == 1:
                    timestepper_flow = list_
                elif process_model == 2:
                    timestepper_tran = list_
            elif block_card.startswith('NEWTON_SOLVER'):
                if process_model == 1:
                    newton_solver_flow = list_
                elif process_model == 2:
                    newton_solver_tran = list_
            elif block_card.startswith('LINEAR_SOLVER'):
                if process_model == 1:
                    linear_solver_flow = list_
                elif process_model == 2:
                    linear_solver_tran = list_
        
        if store_mode:
            list_.append(line)
            
            
    
    append_dict_to_list(timestepper_flow,timestepper_flow_dict)
    ##append transport dic only if other dics aren't empy
    if timestepper_tran or newton_solver_tran or linear_solver_tran:
        append_dict_to_list(timestepper_tran,timestepper_trans_dict)
    append_dict_to_list(newton_solver_flow,newton_flow_dict)
    if timestepper_tran or newton_solver_tran or linear_solver_tran:
        append_dict_to_list(newton_solver_tran,newton_trans_dict)

    f.seek(0)
    f2 = open(filename+'.tmp','w')
#    print('here')
    skip_mode = False
    insert_mode = False
    while True:
        line = f.readline()
        if len(line) == 0:
            break
        w = line.split()
        card = ''
        if len(w) > 0:
            card = w[0].strip().upper()
#        print(card)
        if card.startswith(tupl):
            skip_mode = True
        if not skip_mode:
#            if card not in timestepper_flow_dict and timestepper_trans_dict \
#               and newton_flow_dict and newton_trans_dict:
#                   f2.write(line)
            if card not in newton_flow_dict and card not in timestepper_flow_dict \
              and card not in newton_trans_dict and card not in timestepper_trans_dict:
                   f2.write(line)
        if card.startswith('END') or card.startswith('/'):
            skip_mode = False
        if card == 'SUBSURFACE':
            insert_mode = True

        if insert_mode:
            if timestepper_flow or newton_solver_flow or linear_solver_flow or \
                timestepper_tran or newton_solver_tran or linear_solver_tran:
                    
                f2.write('\n#=========================== numerical methods '
                     '================================\n')
                if timestepper_flow or newton_solver_flow or linear_solver_flow:
                    f2.write('NUMERICAL_METHODS FLOW\n')

                    if timestepper_flow:
                        write_block(f2,'TIMESTEPPER',timestepper_flow)
                    if newton_solver_flow:
                        write_block(f2,'NEWTON_SOLVER',newton_solver_flow)   
                    if linear_solver_flow:
                        write_block(f2,'LINEAR_SOLVER',linear_solver_flow)

                if timestepper_tran or newton_solver_tran or linear_solver_tran:

                    f2.write('NUMERICAL_METHODS TRANSPORT\n')
                    if timestepper_tran:
                        write_block(f2,'TIMESTEPPER',timestepper_tran)
                    if newton_solver_tran:
                        write_block(f2,'NEWTON_SOLVER',newton_solver_tran)
                    if linear_solver_tran:
                        write_block(f2,'LINEAR_SOLVER',linear_solver_tran)

                f2.write('END')
            insert_mode = False
    f.close()
    f2.close()
    
    if replace_file_flag:
        print('File {} has been updated.'.format(filename))
        # using shutil.move adds ^M to end of lines.
        os.remove(filename)
        shutil.copy(filename+'.tmp',filename)
#    os.remove(filename+'.tmp')

replace_file_flag = False


def main():
#    filename = '1d_flux.in'
    suffix = '*.in'
    for root, dirnames, filenames in os.walk('.'):
        for filename in fnmatch.filter(filenames,suffix):
            filename = os.path.join(root,filename)
            print(filename)
            refactor_file(filename,replace_file_flag)


if __name__ == "__main__":
    try:
        suite_status = main()
        print("success")
        sys.exit(suite_status)
    except Exception as error:
        print(str(error))
#        if cmdl_options.backtrace:
#            traceback.print_exc()
#        traceback.print_exc()
        print("failure")
        sys.exit(1)

print('done')
