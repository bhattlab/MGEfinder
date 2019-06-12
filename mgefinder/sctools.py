import warnings
warnings.filterwarnings("ignore")

def left_softclip_length(read):
    return read.cigartuples[0][1]

def right_softclip_length(read):
    return read.cigartuples[-1][1]

def left_softclipped_position(read):
    return read.get_reference_positions()[0]-1

def right_softclipped_position(read):
    return read.get_reference_positions()[-1] + 1

def is_right_softclipped_lenient(read):
    if read.cigartuples[-1][0] == 4:
        return True
    elif read.get_tag('MD')[-1] == '0' and not read.get_tag('MD')[-2].isdigit():
        return True
    else:
        return False

def is_left_softclipped_lenient(read):
    if read.cigartuples[0][0] == 4:
        return True
    elif read.get_tag('MD')[0] == '0':
        return True
    else:
        return False

def is_left_softclipped_lenient_at_site(read, contig, pos):
    if not is_left_softclipped_lenient(read):
        return False
    if left_softclipped_site_lenient(read) == (contig, pos):
        return True
    else:
        return False

def is_right_softclipped_lenient_at_site(read, contig, pos):
    if not is_right_softclipped_lenient(read):
        return False
    if right_softclipped_site_lenient(read) == (contig, pos):
        return True
    else:
        return False

def is_softclipped_lenient_at_site(read, contig, pos):
    if is_right_softclipped_lenient_at_site(read, contig, pos):
        return True
    elif is_left_softclipped_lenient_at_site(read, contig, pos):
        return True
    else:
        return False

def is_left_softclipped_strict(read):
    if read.cigartuples[0][0] == 4:
        return True
    else:
        return False

def is_right_softclipped_strict(read):
    if read.cigartuples[-1][0] == 4:
        return True
    else:
        return False

def get_right_softclip_length(read):
    if is_right_softclipped_lenient(read):
        if read.cigartuples[-1][0] == 4:
            return read.cigartuples[-1][1]
        else:
            return 1
    else:
        return 0

def get_right_softclip_length_strict(read):
    if is_right_softclipped_strict(read):
        return read.cigartuples[-1][1]
    else:
        return 0

def get_left_softclip_length(read):
    if is_left_softclipped_lenient(read):
        if read.cigartuples[0][0] == 4:
            return read.cigartuples[0][1]
        else:
            return 1
    else:
        return 0

def get_left_softclip_length_strict(read):
    if is_left_softclipped_strict(read):
        return read.cigartuples[0][1]
    else:
        return 0

def right_softclipped_site_lenient(read):
    if read.cigartuples[-1][0] == 4:
        return read.reference_name, read.get_reference_positions()[-1] + 1
    elif read.query_sequence[-1] != read.get_reference_sequence()[-1]:
        return read.reference_name, read.get_reference_positions()[-1]

def left_softclipped_site_lenient(read):
    if read.cigartuples[0][0] == 4:
        return read.reference_name, read.get_reference_positions()[0] - 1
    elif read.query_sequence[0] != read.get_reference_sequence()[0]:
        return read.reference_name, read.get_reference_positions()[0]

def right_softclipped_sequence(read):
    if is_right_softclipped_lenient(read):
        if read.cigartuples[-1][0] == 4:
            return read.query_sequence[-read.cigartuples[-1][1]:]
        else:
            return read.query_sequence[-1]
    else:
        return ''

def right_softclipped_sequence_strict(read):
    if is_right_softclipped_strict(read):
        return read.query_sequence[-read.cigartuples[-1][1]:]
    else:
        return ''

def right_softclip_qualities(read):
    if is_right_softclipped_lenient(read):
        if read.cigartuples[-1][0] == 4:
            return list(read.query_qualities)[-read.cigartuples[-1][1]:]
        else:
            return [list(read.query_qualities)[-1]]
    else:
        return []

def left_softclipped_sequence(read):
    if is_left_softclipped_lenient(read):
        if read.cigartuples[0][0] == 4:
            return read.query_sequence[:read.cigartuples[0][1]]
        else:
            return read.query_sequence[0]
    else:
        return ''

def left_softclipped_sequence_strict(read):
    if is_left_softclipped_strict(read):
        return read.query_sequence[:read.cigartuples[0][1]]
    else:
        return ''

def left_softclip_qualities(read):
    if is_left_softclipped_lenient(read):
        if read.cigartuples[0][0] == 4:
            return list(read.query_qualities)[:read.cigartuples[0][1]]
        else:
            return [list(read.query_qualities)[0]]
    else:
        return []

def is_double_softclipped_lenient(read):

    if is_left_softclipped_lenient(read) and is_right_softclipped_lenient(read):
        return True
    else:
        return False

def read_meets_min_alignment_inner_length(read, min_alignment_inner_length):

    if not is_double_softclipped_lenient(read):
        return True

    right_length = len(right_softclipped_sequence(read))
    left_length = len(left_softclipped_sequence(read))

    total_length = len(read.query_sequence)
    aligned_length = total_length - right_length - left_length

    if aligned_length >= min_alignment_inner_length:
        return True
    else:
        return False

def right_softclip_proportion(read):
    return len(right_softclipped_sequence_strict(read)) / len(read.query_sequence)

def left_softclip_proportion(read):
    return len(left_softclipped_sequence_strict(read)) / len(read.query_sequence)

def left_softclip_reference_start(read):
    return read.reference_start - len(left_softclipped_sequence_strict(read))

def right_softclip_reference_end(read):
    return read.reference_end + len(right_softclipped_sequence_strict(read))

