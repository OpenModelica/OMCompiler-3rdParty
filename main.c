/* Test program for trying to link dgesv (used to figure out which
 * files are actually needed to link against only dgesv).
 */
int dgesv_();
int main()
{
  return dgesv_();
}
