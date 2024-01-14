
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

################# Алгоритм Нидлмана Вунша #################

class Lists():
    
    instances = []
    def __init__(self, dct, elem_ind):
        self.dct = dct
        self.elem_ind = elem_ind
        self.__class__.instances.append(self)

        
def nid_wun(seq_1, seq_2):
    
    match = 1
    mismatch = -1
    gap = -1
    
    nw_matrix = []
    for i in range(len(seq_1) + 1):
        nw_matrix.append([])
    

    # инициализация таблицы
    
    nw_matrix[0].append(0)
    
    # рассчитываем значения первого ряда
    for i in range(len(seq_1)):
        nw_matrix[0].append(-i + -1)

    
    # рассчитываем значения первого столбца
    for i, j in enumerate(nw_matrix[1:]):
        j.append(-i + -1)
    

    # заполняем таблицу
    backarrows = {}
    
    for m, symb_1 in enumerate(seq_1):
        for n, symb_2 in enumerate(seq_2):

    
            if symb_1 == symb_2:
                max_elem = nw_matrix[m][n] + match
                max_elem_ind = [[m, n, 1]]
                
            else:
                left_elem = nw_matrix[m + 1][n] + gap
                topleft_elem = nw_matrix[m][n] + mismatch
                top_elem = nw_matrix[m][n + 1] + gap 

                max_left_ind = [m + 1, n, 0]
                max_topleft_ind = [m, n, 1]
                max_top_ind = [m, n + 1, 2]
                
                surr = [left_elem, topleft_elem, top_elem]
                surr_ind = [max_left_ind, max_topleft_ind, max_top_ind]
                
                max_elem = max(surr)
                max_list = [i for i, j in enumerate(surr) if j == max_elem]
                max_elem_ind =  [j for i, j in enumerate(surr_ind) if i in max_list]
                
                
            backarrows[m + 1, n + 1] = max_elem_ind
            
            nw_matrix[m + 1].append(max_elem)

    
    # возвращаемся и создаем варианты выравнивания
    elem_ind = (len(seq_1),len(seq_2))
    list_1 = Lists({}, elem_ind) # Первая ветвь выравнивания
   
    
    def backgoing(branch_obj, backarrows, count):
        
        import copy
        
        while branch_obj.elem_ind != (0, 0):
            
            if len(backarrows.get(branch_obj.elem_ind)) == 1: # если только одна стрелка
                if backarrows.get(branch_obj.elem_ind)[0][2] == 1: # если стрелка по диагонали
                    
                    if seq_1[branch_obj.elem_ind[0] - 1] != seq_2[branch_obj.elem_ind[1] - 1]: # выясняем матч или мисматч
                        branch_obj.dct[count] = [seq_1[branch_obj.elem_ind[0] - 1], seq_2[branch_obj.elem_ind[1] - 1], mismatch]
                        
                    else:
                        branch_obj.dct[count] = [seq_1[branch_obj.elem_ind[0] - 1], seq_2[branch_obj.elem_ind[1] - 1], match]
                        
                    branch_obj.elem_ind = (backarrows.get(branch_obj.elem_ind)[0][0], backarrows.get(branch_obj.elem_ind)[0][1]) # замена на следующий элемент
                
                
                    
                elif backarrows.get(branch_obj.elem_ind)[0][2] == 0: # если стрелка влево
                    branch_obj.dct[count] = [seq_2[branch_obj.elem_ind[1] - 1], '-', gap]
                    
                    branch_obj.elem_ind = (backarrows.get(branch_obj.elem_ind)[0][0], backarrows.get(branch_obj.elem_ind)[0][1]) # замена на новый элемент
                    
                
                elif backarrows.get(branch_obj.elem_ind)[0][2] == 2: # если стрелка вверх
                    branch_obj.dct[count] = [seq_1[branch_obj.elem_ind[0] - 1], '-', gap]
                    
                    branch_obj.elem_ind = (backarrows.get(branch_obj.elem_ind)[0][0], backarrows.get(branch_obj.elem_ind)[0][1]) # замена на новый элемент
            
            
                elif backarrows.get(branch_obj.elem_ind)[0][2] == 3:
                    branch_obj.dct[count] = [seq_1[branch_obj.elem_ind[0] - 1], seq_2[branch_obj.elem_ind[1] - 1], mismatch]
                    
                    
                    branch_obj.elem_ind = (backarrows.get(branch_obj.elem_ind)[0][0], backarrows.get(branch_obj.elem_ind)[0][1])
                    backarrows[branch_obj.elem_ind][0][2] = 2
                
                
                elif backarrows.get(branch_obj.elem_ind)[0][2] == 4:
                    branch_obj.dct[count] = [seq_1[branch_obj.elem_ind[0] - 1], seq_2[branch_obj.elem_ind[1] - 1], mismatch]
                    
                    branch_obj.elem_ind = (backarrows.get(branch_obj.elem_ind)[0][0], backarrows.get(branch_obj.elem_ind)[0][1])
                    backarrows[branch_obj.elem_ind][0][2] = 0
                
                
            elif len(backarrows.get(branch_obj.elem_ind)) > 1: # если больше одной стрелки
                arrows = backarrows.get(branch_obj.elem_ind) # все возможные стрелки
                
                new_dicts = []
                for index, i in enumerate(arrows): 
                    new_backarrows = copy.deepcopy(backarrows)
                    
                    for ind, j in enumerate(arrows):
                        if i[2] != j[2]:
                            del new_backarrows[branch_obj.elem_ind][ind]
                    
                    if new_backarrows[branch_obj.elem_ind][0][2] == 2:
                        new_backarrows[branch_obj.elem_ind][0][2] = 3
                    if new_backarrows[branch_obj.elem_ind][0][2] == 0:
                        new_backarrows[branch_obj.elem_ind][0][2] = 4
                        
                    new_dicts.append(new_backarrows)
                    
                for i in new_dicts:
                    new_branch = Lists(copy.deepcopy(branch_obj.dct), branch_obj.elem_ind)
                    backgoing(new_branch, i, count)
                break
            
        
            count += 1
            
    backgoing(list_1, backarrows, 0)

    
    for i in Lists.instances:
        sequence_1 = []
        sequence_2 = []
        results = []
        results_2 = []
        score = 0
        for j in i.dct.values():
            sequence_1.append(j[0])
            sequence_2.append(j[1])
            results.append(str(j[2]))
        for i in results:
            if i != '1':
                results_2.append(' ')
                score -= 1
            else:
                results_2.append('|')
                score += 1
        print(''.join(sequence_1), end='\n')
        print(''.join(results_2), end='\n')
        print(''.join(sequence_2), end='\n')
        print(f'Score: {score}', end='\n\n')
    

nid_wun('GATTACA', 'GCATGCG')
x = pairwise2.align.globalms('GCATGCG', 'GATTACA', 1, -1, -1, 0)
print(format_alignment(*x[0]))
