/* do we need to process this record? 1=yes 0=no -1=break*/
int process(float time, float start_time, float final_time) /* includefile */
{
  if (start_time == 0.0 && final_time == 0.0) {
    return 1;
  } else if (start_time == 0.0 && time <= final_time) {
    return 1;
  } else if (final_time == 0.0 && time >= start_time) {
    return 1;
  } else if (time >= start_time && time <= final_time) {
    return 1;
  } else if (final_time == 0.0 && time < start_time) {
    return 0;
  } else if (time >= final_time) {
    return -1;
  } else {
    return 0;
  }
}
